"""Le modèle tel que la validation le lit, et les trois références de `docs/S21_CRITERES.md`.

- **Proximité** : distance médiane entre billes, sur les repliements et les paires de copies.
- **Contact** : centres à moins de 1,5 × (r_i + r_j), la définition de la semaine 7.
- **N1** (distance seule) : la moyenne du modèle à la même séparation, sans la 3D.
- **N2** (permutations) : identités des billes permutées à l'intérieur de chaque chromosome.
- **Témoin positif** : une moitié des repliements contre l'autre.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from ..ensemble import store

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
CUTOFF = 1.5            # contact : centres à moins de 1,5 × (r_i + r_j) — semaine 7
PERMUTATIONS = 1000
SEED = 21


@dataclass
class Model:
    """Le noyau et ses repliements, indexés pour la validation."""

    coords: np.ndarray        # (S, n, 3) µm, float32
    radius: np.ndarray        # (n,) µm
    copy_id: np.ndarray       # (n,)
    start: np.ndarray         # (n,) hg38, 0-based
    end: np.ndarray
    labels: tuple[str, ...]   # étiquette de chaque copie : « chr1:a », « chrX:Xi »…
    nuclear_radius: float
    lad_source: str
    lad: np.ndarray           # la piste LAD (synthétique) que le noyau a reçue

    @classmethod
    def load(cls, path: Path | None = None) -> "Model":
        e = store.read(path or DATA / "ensemble" / "gm12878.zarr")
        b = e.beads
        return cls(coords=e.coords.astype(np.float32), radius=b.radius, copy_id=b.copy_id,
                   start=b.start, end=b.end, labels=b.labels, nuclear_radius=b.nuclear_radius,
                   lad_source=b.lad_source, lad=b.lad)

    def subset(self, structures: np.ndarray) -> "Model":
        return Model(self.coords[structures], self.radius, self.copy_id, self.start, self.end,
                     self.labels, self.nuclear_radius, self.lad_source, self.lad)

    def chrom_of_copy(self, c: int) -> str:
        return self.labels[c].split(":", 1)[0]

    def copies(self, chrom: str) -> list[int]:
        return [c for c in range(len(self.labels)) if self.chrom_of_copy(c) == chrom]

    def beads_at(self, chrom: str, pos: np.ndarray) -> np.ndarray:
        """(copies, k) : la bille de chaque copie qui contient chaque position ; −1 sinon."""
        out = []
        for c in self.copies(chrom):
            ids = np.flatnonzero(self.copy_id == c)
            s, e = self.start[ids], self.end[ids]
            k = np.searchsorted(s, pos, side="right") - 1
            ok = (k >= 0) & (k < len(ids))
            k = np.clip(k, 0, len(ids) - 1)
            ok &= pos < e[k]
            out.append(np.where(ok, ids[k], -1))
        return np.array(out, dtype=np.int64)

    def radial(self) -> np.ndarray:
        """(n,) position radiale moyenne sur les repliements, en fraction du rayon du noyau."""
        return (np.linalg.norm(self.coords, axis=2) / self.nuclear_radius).mean(axis=0)


def pair_distances(model: Model, idx: np.ndarray, *, cis_only: bool, cols: np.ndarray | None = None) -> np.ndarray:
    """(échantillons, k, k) distances entre régions, pour chaque repliement et paire de copies.

    `idx` vient de `beads_at` : (copies, k). Une région non placée sur une copie donne NaN.
    `cols`, de même forme, donne les colonnes (sinon les mêmes régions) : (échantillons, k, k′).
    `cis_only` : seulement les paires dans la même copie (traçage, copie par copie) ; sinon
    toutes les paires de copies (données non phasées : SPRITE).
    """
    X = model.coords
    other = idx if cols is None else cols
    combos = [(a, a) for a in range(len(idx))] if cis_only else [(a, b) for a in range(len(idx)) for b in range(len(other))]
    out = []
    for a, b in combos:
        ia, ib = idx[a], other[b]
        pa = X[:, np.maximum(ia, 0)]
        pb = X[:, np.maximum(ib, 0)]
        d = np.linalg.norm(pa[:, :, None, :] - pb[:, None, :, :], axis=3)
        d[:, ia < 0, :] = np.nan
        d[:, :, ib < 0] = np.nan
        out.append(d)
    return np.concatenate(out, axis=0)


def pair_radii(model: Model, idx: np.ndarray) -> np.ndarray:
    """(k, k) somme des rayons de deux régions (copie a), pour le seuil de contact."""
    r = model.radius[np.maximum(idx[0], 0)]
    return r[:, None] + r[None, :]


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    ra, rb = rankdata(a[ok]), rankdata(b[ok])
    return float(np.corrcoef(ra, rb)[0, 1])


def over_expected(m: np.ndarray, offsets: np.ndarray) -> np.ndarray:
    """Observé / attendu : chaque valeur divisée par la moyenne de sa diagonale (`offsets`)."""
    out = np.full(m.shape, np.nan)
    for k in np.unique(offsets[np.isfinite(m)]):
        sel = (offsets == k) & np.isfinite(m)
        mu = m[sel].mean()
        if mu != 0:
            out[sel] = m[sel] / mu
    return out


def expected(m: np.ndarray, offsets: np.ndarray) -> np.ndarray:
    """N1 : la moyenne de la diagonale à la place de chaque valeur — la distance seule."""
    out = np.full(m.shape, np.nan)
    for k in np.unique(offsets[np.isfinite(m)]):
        sel = (offsets == k) & np.isfinite(m)
        out[sel] = m[sel].mean()
    return out


@dataclass
class Result:
    """Un critère : la valeur mesurée, sa référence, son verdict."""

    id: str
    data: str
    measure: str
    value: float
    threshold: str
    verdict: str
    null_p95: float | None = None
    reference: float | None = None
    extra: dict | None = None


def verdict(value: float, threshold: float, *, sign: int = 1, p95: float | None = None) -> str:
    """« accord » au-delà du seuil (et des permutations), « partiel » du bon signe et au-delà des
    permutations mais sous le seuil, « échec » sinon — `docs/S21_CRITERES.md`."""
    if not np.isfinite(value):
        return "échec"
    v = sign * value
    beyond = p95 is None or v > p95
    if v >= sign * threshold and beyond:
        return "accord"
    if v > 0 and beyond and p95 is not None:
        return "partiel"
    return "échec"


def null_p95(stat, permute, *, n: int = PERMUTATIONS, seed: int = SEED, sign: int = 1) -> tuple[float, float]:
    """95ᵉ centile (et moyenne) de `sign × stat(permute(rng))` sur `n` permutations."""
    rng = np.random.default_rng(seed)
    vals = np.array([sign * stat(permute(rng)) for _ in range(n)])
    vals = vals[np.isfinite(vals)]
    return float(np.percentile(vals, 95)), float(vals.mean())


def halves(model: Model) -> tuple[Model, Model]:
    s = np.arange(model.coords.shape[0])
    return model.subset(s[0::2]), model.subset(s[1::2])
