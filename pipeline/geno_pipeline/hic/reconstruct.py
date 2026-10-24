"""Contact → distance → coordonnées 3D.

La conversion la plus critique du projet, et la plus facile à faire passer pour
acquise. Toute la chaîne repose sur une loi de puissance :

    d(i, j) ∝ f(i, j)^(-alpha)

`alpha` n'est pas une constante de la nature. C'est un paramètre, souvent repris
à 1/3 de proche en proche sans être revérifié. Si le modèle direct décroît en
d^(-gamma), l'exposant correct est 1/gamma — et rien ne garantit que gamma vaille
3 dans le jeu qu'on a sous la main.

Trois pièges, traités explicitement ici :

1. **Les contacts nuls ne sont pas des distances infinies.** Deux billes sans
   contact observé sont loin, pas infiniment loin. ShRec3D (Lesne 2014) complète
   par les plus courts chemins dans le graphe des contacts : la distance devient
   géodésique, ce qui est à la fois fini et cohérent avec l'inégalité triangulaire.

2. **Une matrice de distances ne détermine la structure qu'à une isométrie près.**
   Rotation, translation, et surtout **réflexion** : la chiralité est
   irrécupérable depuis des distances seules. Toute comparaison à une vérité doit
   donc passer par un alignement de Procruste qui autorise la réflexion. Sans ça,
   on mesure l'orientation du repère, pas la qualité de la reconstruction.

3. **Le RMSD brut dépend de l'échelle**, elle-même arbitraire. On le normalise
   par le rayon de giration de la vérité.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def balance_dense(counts: np.ndarray, iters: int = 60, tol: float = 1e-5) -> np.ndarray:
    """Équilibrage ICE sur une matrice dense : marginales rendues uniformes.

    Même principe que `cooler.balance_cooler`, en quelques lignes, pour les
    matrices en mémoire du banc d'essai. Sans ça, un biais de couverture se lit
    comme une proximité spatiale.
    """
    m = counts.astype(np.float64).copy()
    n = len(m)
    weights = np.ones(n)
    for _ in range(iters):
        marg = m.sum(axis=1)
        nz = marg > 0
        if not nz.any():
            break
        scale = np.ones(n)
        scale[nz] = marg[nz] / marg[nz].mean()
        weights[nz] /= scale[nz]
        m /= np.outer(np.where(nz, scale, 1.0), np.where(nz, scale, 1.0))
        if np.nanstd(marg[nz] / marg[nz].mean()) < tol:
            break
    return m


def to_distance(counts: np.ndarray, alpha: float) -> np.ndarray:
    """d = f^(-alpha). Les contacts nuls deviennent `inf`, à compléter ensuite."""
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where(counts > 0, counts.astype(np.float64) ** (-alpha), np.inf)
    np.fill_diagonal(d, 0.0)
    return d


def complete_shortest_path(dist: np.ndarray) -> np.ndarray:
    """Remplace les distances manquantes par la géodésique du graphe.

    C'est l'apport de ShRec3D. Une paire sans contact observé n'est pas à
    l'infini : elle est au mieux à la somme des sauts qui la relient. Le résultat
    respecte l'inégalité triangulaire, ce qu'une matrice trouée ne fait pas — et
    le MDS classique, lui, la suppose.
    """
    from scipy.sparse.csgraph import shortest_path

    finite = np.where(np.isfinite(dist), dist, 0.0)
    geo = shortest_path(finite, method="D", directed=False)
    if not np.isfinite(geo).all():                  # graphe déconnecté
        geo = np.where(np.isfinite(geo), geo, np.nanmax(geo[np.isfinite(geo)]) * 2)
    np.fill_diagonal(geo, 0.0)
    return geo


def classical_mds(dist: np.ndarray, dim: int = 3) -> np.ndarray:
    """MDS métrique (Torgerson) : matrice de distances → coordonnées.

    Double centrage pour obtenir la matrice de Gram, puis ses vecteurs propres
    dominants. Les valeurs propres négatives signalent une matrice non
    euclidienne — inévitable après complétion — et sont écrêtées à zéro.
    """
    n = len(dist)
    d2 = dist**2
    j = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * j @ d2 @ j
    gram = (gram + gram.T) / 2                      # symétrie exacte
    vals, vecs = np.linalg.eigh(gram)
    order = np.argsort(vals)[::-1][:dim]
    lam = np.clip(vals[order], 0, None)
    return vecs[:, order] * np.sqrt(lam)


def procrustes(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    """Aligne `x` sur `y` — translation, rotation, **réflexion**, échelle.

    La réflexion est autorisée délibérément : une matrice de distances ne porte
    aucune information de chiralité, donc une reconstruction miroir est une
    reconstruction correcte. Interdire la réflexion reviendrait à compter comme
    erreur la moitié des solutions valides.
    """
    xc = x - x.mean(axis=0)
    yc = y - y.mean(axis=0)
    u, s, vt = np.linalg.svd(xc.T @ yc)
    rot = u @ vt                                    # pas de correction de det
    denom = (xc**2).sum()
    scale = s.sum() / denom if denom > 0 else 1.0
    aligned = scale * xc @ rot
    rmsd = float(np.sqrt(((aligned - yc) ** 2).sum(axis=1).mean()))
    return aligned, rmsd


@dataclass(frozen=True)
class Fidelity:
    """Qualité d'une reconstruction face à une structure connue."""

    alpha: float
    nrmsd: float           # RMSD après alignement, en rayons de giration
    dist_rho: float        # Spearman sur les distances deux à deux
    radial_r: float        # Pearson sur la position radiale

    def __str__(self) -> str:
        return (
            f"alpha {self.alpha:.3f}   nRMSD {self.nrmsd:.4f}   "
            f"rho(distances) {self.dist_rho:+.4f}   r(radial) {self.radial_r:+.4f}"
        )


def _gyration(coords: np.ndarray) -> float:
    return float(np.sqrt(((coords - coords.mean(axis=0)) ** 2).sum(axis=1).mean()))


def _radial(coords: np.ndarray) -> np.ndarray:
    return np.linalg.norm(coords - coords.mean(axis=0), axis=1)


def reconstruct(counts: np.ndarray, alpha: float, dim: int = 3) -> np.ndarray:
    """Chaîne complète : comptages → distances → complétion → coordonnées."""
    return classical_mds(complete_shortest_path(to_distance(counts, alpha)), dim)


def fidelity(recon: np.ndarray, truth: np.ndarray, alpha: float) -> Fidelity:
    from scipy.stats import pearsonr, spearmanr

    _aligned, rmsd = procrustes(recon, truth)
    gyr = _gyration(truth)

    iu = np.triu_indices(len(truth), 1)
    d_true = np.linalg.norm(truth[:, None] - truth[None, :], axis=2)[iu]
    d_rec = np.linalg.norm(recon[:, None] - recon[None, :], axis=2)[iu]

    return Fidelity(
        alpha=alpha,
        nrmsd=rmsd / gyr if gyr > 0 else float("inf"),
        dist_rho=float(spearmanr(d_true, d_rec).statistic),
        radial_r=float(pearsonr(_radial(truth), _radial(recon)).statistic),
    )


def sweep(counts: np.ndarray, truth: np.ndarray, alphas) -> list[Fidelity]:
    """Balaie alpha et rend la fidélité à chaque valeur."""
    return [fidelity(reconstruct(counts, a), truth, a) for a in alphas]


def best(results: list[Fidelity]) -> Fidelity:
    """Le meilleur alpha au sens du nRMSD — la métrique qui juge la géométrie,
    pas seulement l'ordre des distances."""
    return min(results, key=lambda f: f.nrmsd)
