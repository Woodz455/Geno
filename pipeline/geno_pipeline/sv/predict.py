"""Ce que l'extrusion 1D seule dit des contacts — assez vite pour un navigateur.

Une boucle tenue par une cohésine rapproche ses deux pieds : sur la chaîne, ils sont à un pas
l'un de l'autre. La **distance effective** entre deux monomères est donc le plus court chemin le
long de la chaîne, chaque boucle tenue servant de raccourci — un graphe de `n` nœuds et de
`n − 1 + boucles` arêtes de longueur 1, qu'un parcours en largeur résout. Les boucles ne se
croisent pas (deux jambes ne se dépassent pas), mais le parcours n'en a pas besoin.

La probabilité de contact s'en déduit par une seule loi, `P = min(1, (d / d0)^−α + c)`, moyennée
sur les instantanés : une décroissance, et un plancher `c` — le mélange qu'impose le confinement,
que le 3D montre au-delà d'un mégabase (~0,1 % de contacts quelle que soit la distance). `α` et `d0` ne sont pas choisis : ils sont **calés sur le 3D du type sauvage**
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


def law(d: np.ndarray, alpha: float, d0: float, c: float = 0.0) -> np.ndarray:
    with np.errstate(divide="ignore"):
        p = np.where(d > 0, (np.maximum(d, 1e-9) / d0) ** -alpha + c, 1.0)
    return np.minimum(1.0, p)


def rows(n: int, legs: np.ndarray, sources: np.ndarray, alpha: float, d0: float, c: float = 0.0) -> np.ndarray:
    """(len(sources), n) : probabilité de contact prédite, moyennée sur les instantanés."""
    acc = np.zeros((len(sources), n))
    for lg in legs:
        acc += law(shortest_path(graph(n, lg), method="D", directed=False, unweighted=True, indices=sources),
                   alpha, d0, c)
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
    c: float              # plancher : le mélange dû au confinement
    r_log: float          # Pearson des log, casiers de 20 kb (10 monomères), séparation ≥ 3
    pairs: int


def binned(rows_: np.ndarray, k: int = 10) -> np.ndarray:
    """Lignes sommées par casiers de `k` monomères (20 kb) : à 2 kb, une paire a rarement plus
    d'un contact sur 1 200 instantanés, et comparer des log de zéros ne mesure rien."""
    m = rows_.shape[1] // k * k
    return rows_[:, :m].reshape(rows_.shape[0], -1, k).mean(axis=2)


def calibrate(coords: np.ndarray, legs: np.ndarray, *, every: int = 25, frames: int = 300,
              cutoff: float = 1.5, floor: bool = True, max_sep: int | None = None) -> Fit:
    """Cale `α, d0, c` sur une exécution 3D, par la **déviance de Poisson** sur toutes les paires,
    zéros compris : le contact observé d'une paire est un compte sur les instantanés, et la
    première version — moindres carrés sur le log des seules paires non nulles (20 %) — sélectionnait
    les paires chanceuses et convergeait vers une loi dégénérée (α = 0,23, d0 → 0, r = 0,61).

    Contacts observés : tous les instantanés ; distances effectives : `frames` d'entre eux,
    réguliers (la moyenne de la loi n'a pas besoin de plus). Grille grossière, puis Nelder–Mead
    sur `(α, log d0, log c)`.

    `floor=False, max_sep` (révision) : sans plancher (`c = 0`), calé sur les seules paires à
    moins de `max_sep` monomères — là où la sphère de confinement ne décide pas encore."""
    from scipy.optimize import minimize

    n = coords.shape[1]
    src = np.arange(0, n, every)
    obs = observed(coords, src, cutoff)
    pick = np.linspace(0, len(coords) - 1, min(frames, len(coords))).astype(int)
    d = distances(n, legs[pick], src)                                # (f, s, n)
    sep = np.abs(src[:, None] - np.arange(n)[None, :])
    ok = (sep >= 3) & (sep < (max_sep or n))
    F = len(coords)
    o = obs[ok]

    def full(x):
        return x if floor else np.r_[x[:2], -np.inf]

    def dev(x) -> float:
        x = full(x)
        p = np.clip(law(d, x[0], float(np.exp(x[1])), float(np.exp(x[2]))).mean(axis=0)[ok], 1e-12, 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(o > 0, o * np.log(o / p), 0.0)
        return float(2 * F * np.sum(t - (o - p)) / o.size)

    grid = [(a, np.log(d0), np.log(c)) for a in (1.0, 1.5, 2.0) for d0 in (0.5, 1.0, 2.0) for c in (3e-4, 1e-3, 3e-3)]
    if not floor:
        grid = sorted({g[:2] for g in grid})
    x0 = min(grid, key=dev)
    r = minimize(dev, np.array(x0), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-7, "maxiter": 400})
    x = full(r.x)
    a, d0, c = float(x[0]), float(np.exp(x[1])), float(np.exp(x[2]))
    p = law(d, a, d0, c).mean(axis=0)
    bo, bp = binned(np.where(ok, obs, np.nan)), binned(np.where(ok, p, np.nan))
    good = np.isfinite(bo) & (bo > 0)
    rr = float(np.corrcoef(np.log(bp[good]), np.log(bo[good]))[0, 1])
    return Fit(round(a, 4), round(d0, 4), float(f"{c:.4g}"), rr, int(ok.sum()))


def contact(p_rows: np.ndarray, src: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """Contact moyen entre deux groupes de monomères, depuis les lignes de sources `src`
    (qui doivent contenir `a`)."""
    at = {int(s): k for k, s in enumerate(src)}
    return float(np.mean([p_rows[at[int(i)], b].mean() for i in a]))


__all__ = ["Fit", "binned", "calibrate", "contact", "contact_map", "distances", "graph", "law", "observed", "rows"]
