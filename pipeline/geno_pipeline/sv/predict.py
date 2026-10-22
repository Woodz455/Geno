"""Ce que l'extrusion 1D seule dit des contacts — assez vite pour un navigateur.

Une boucle tenue par une cohésine rapproche ses deux pieds : sur la chaîne, ils sont à un pas
l'un de l'autre. La **distance effective** entre deux monomères est donc le plus court chemin le
long de la chaîne, chaque boucle tenue servant de raccourci — un graphe de `n` nœuds et de
`n − 1 + boucles` arêtes de longueur 1, qu'un parcours en largeur résout. Les boucles ne se
croisent pas (deux jambes ne se dépassent pas), mais le parcours n'en a pas besoin.

La probabilité de contact s'en déduit par une seule loi, `P = min(1, (d / d0)^−α)`, moyennée sur
les instantanés. `α` et `d0` ne sont pas choisis : ils sont **calés sur le 3D du type sauvage**
(OpenMM, mêmes jambes), puis la loi est **éprouvée** sur le 3D des allèles mutants, qu'elle n'a
pas vus. Ce qu'elle ignore — volume exclu, confinement, enchevêtrement — est ce que ce test
mesure.

Le port TypeScript (`apps/variant/extrusion.ts`) applique les mêmes règles et la même loi ; la
mesure le confronte à celui-ci.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from ..fine.observe import contact_map


def graph(n: int, legs: np.ndarray) -> sparse.csr_matrix:
    """La chaîne et les boucles tenues d'un instantané (`legs` : (n_lef, 2), −1 = décroché)."""
    held = legs[(legs[:, 0] >= 0) & (legs[:, 1] >= 0) & (legs[:, 0] != legs[:, 1])]
    i = np.r_[np.arange(n - 1), held[:, 0]]
    j = np.r_[np.arange(1, n), held[:, 1]]
    g = sparse.coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)).tocsr()
    return g


def distances(n: int, legs: np.ndarray, sources: np.ndarray) -> np.ndarray:
    """(frames, len(sources), n) : distance effective, en monomères."""
    out = np.empty((len(legs), len(sources), n), np.float32)
    for f, lg in enumerate(legs):
        out[f] = shortest_path(graph(n, lg), method="D", directed=False, unweighted=True, indices=sources)
    return out


def law(d: np.ndarray, alpha: float, d0: float) -> np.ndarray:
    with np.errstate(divide="ignore"):
        p = np.where(d > 0, (np.maximum(d, 1e-9) / d0) ** -alpha, 1.0)
    return np.minimum(1.0, p)


def rows(n: int, legs: np.ndarray, sources: np.ndarray, alpha: float, d0: float) -> np.ndarray:
    """(len(sources), n) : probabilité de contact prédite, moyennée sur les instantanés."""
    acc = np.zeros((len(sources), n))
    for lg in legs:
        acc += law(shortest_path(graph(n, lg), method="D", directed=False, unweighted=True, indices=sources),
                   alpha, d0)
    return acc / max(len(legs), 1)


def observed(coords: np.ndarray, sources: np.ndarray, cutoff: float = 1.5) -> np.ndarray:
    """(len(sources), n) : fréquence de contact 3D (seuil de la semaine 7, en diamètres)."""
    out = np.zeros((len(sources), coords.shape[1]))
    for x in coords:
        d = np.linalg.norm(x[sources][:, None, :] - x[None, :, :], axis=2)
        out += d < cutoff
    out /= len(coords)
    out[np.arange(len(sources)), sources] = 1.0
    return out


@dataclass(frozen=True)
class Fit:
    alpha: float
    d0: float
    r_log: float          # Pearson des log, paires de séparation ≥ 3 monomères, contacts > 0
    pairs: int


def calibrate(coords: np.ndarray, legs: np.ndarray, *, every: int = 20, frames: int = 120,
              cutoff: float = 1.5) -> Fit:
    """Cale `α, d0` sur une exécution 3D : un monomère source tous les `every`, contre toute la
    chaîne, sur `frames` instantanés tirés régulièrement. Grille grossière, puis Nelder–Mead sur
    `(α, log d0)` depuis son meilleur point — une grille fine coûtait des dizaines de minutes."""
    from scipy.optimize import minimize

    n = coords.shape[1]
    pick = np.linspace(0, len(coords) - 1, min(frames, len(coords))).astype(int)
    src = np.arange(0, n, every)
    obs = observed(coords[pick], src, cutoff)
    d = distances(n, legs[pick], src)                                # (f, s, n)
    sep = np.abs(src[:, None] - np.arange(n)[None, :])
    ok = (sep >= 3) & (obs > 0)
    lobs = np.log(obs[ok])

    def loss(x) -> float:
        p = law(d, x[0], float(np.exp(x[1]))).mean(axis=0)
        return float(np.mean((np.log(p[ok]) - lobs) ** 2))

    grid = [(a, np.log(d0)) for a in np.arange(1.0, 2.51, 0.3) for d0 in (0.5, 1.0, 2.0, 4.0)]
    x0 = min(grid, key=loss)
    r = minimize(loss, np.array(x0), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-6, "maxiter": 200})
    a, d0 = float(r.x[0]), float(np.exp(r.x[1]))
    p = law(d, a, d0).mean(axis=0)
    rr = float(np.corrcoef(np.log(p[ok]), lobs)[0, 1])
    return Fit(round(a, 3), round(d0, 4), rr, int(ok.sum()))


def contact(p_rows: np.ndarray, src: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """Contact moyen entre deux groupes de monomères, depuis les lignes de sources `src`
    (qui doivent contenir `a`)."""
    at = {int(s): k for k, s in enumerate(src)}
    return float(np.mean([p_rows[at[int(i)], b].mean() for i in a]))


__all__ = ["Fit", "calibrate", "contact", "contact_map", "distances", "graph", "law", "observed", "rows"]
