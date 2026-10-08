"""Combien pèse un `.g3d` quand le nombre de billes passe de six mille à six cent mille.

Le budget de `ARCHITECTURE.md` § 4 était de l'arithmétique sur des float32. Ce banc le
remplace par des octets **écrits** : on construit le niveau, on l'écrit, on compte.

Deux sortes de structures, et il faut savoir laquelle on regarde :

- **jusqu'à 100 kb par bille**, de vrais noyaux de la semaine 6, construits par le même
  code, mis en cache parce qu'un noyau à 100 kb coûte plusieurs minutes ;
- **en dessous**, un **raffinement synthétique** du noyau à 100 kb : chaque bille est
  découpée en filles le long d'un pont brownien. Ce n'est pas une structure modélisée,
  c'est une chaîne de la bonne longueur et de la bonne rugosité locale. Elle suffit à
  mesurer des octets ; elle ne dit rien de la biologie, et la table le marque.

La variabilité de chaque bille est celle, mesurée en semaine 7, de la bille de 750 kb
qui la contient. La règle de précision est donc la même que dans le vrai fichier.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..nucleus.beads import capacity_at
from ..nucleus.build import build as build_nucleus
from ..fine.region import genome_bp
from .g3d import Level, write

ENSEMBLE = Path(__file__).resolve().parents[2] / "data" / "ensemble" / "gm12878.zarr"
REAL_DOWN_TO = 100_000


def _nucleus(bp: int, cache: Path) -> dict:
    path = cache / f"noyau-{bp}.npz"
    if path.exists():
        z = np.load(path)
        return {k: z[k] for k in z.files}
    nu = build_nucleus(bp_per_bead=bp, seed=0)
    out = {
        "x_nm": nu.x * 1000.0,
        "copy": nu.beads.copy_id,
        "start": nu.beads.start,
        "end": nu.beads.end,
        "radius_nm": nu.beads.radius * 1000.0,
    }
    cache.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **out)
    return out


def _refine(parent: dict, k: int, r_nm: float, seed: int = 0) -> dict:
    """Découpe chaque bille en `k` filles le long d'un pont brownien.

    Le pont part d'un centre de bille et arrive au suivant : la forme à grande échelle
    est celle du noyau réel, seule la rugosité sous la bille est tirée. Le pas des
    filles vaut leur diamètre, ce que la loi de volume de la semaine 6 impose.
    """
    rng = np.random.default_rng(seed)
    xs, copies, starts, ends = [], [], [], []
    copy = parent["copy"]
    for c in np.unique(copy):
        idx = np.flatnonzero(copy == c)
        p = parent["x_nm"][idx]
        nxt = np.r_[p[1:], p[-1:]]
        t = (np.arange(k) / k)[None, :, None]
        base = p[:, None, :] + t * (nxt - p)[:, None, :]            # (m, k, 3)
        walk = np.cumsum(rng.normal(0, 2.0 * r_nm / np.sqrt(3), (len(p), k + 1, 3)), axis=1)
        tt = np.linspace(0, 1, k + 1)[None, :, None]
        bridge = (walk - tt * walk[:, -1:, :])[:, :k, :]          # nul aux deux bouts
        xs.append((base + bridge).reshape(-1, 3))
        s0 = parent["start"][idx]
        span = (parent["end"][idx] - s0) // k
        st = (s0[:, None] + span[:, None] * np.arange(k)[None, :]).reshape(-1)
        starts.append(st)
        ends.append(np.r_[st[1:], parent["end"][idx][-1]])
        copies.append(np.full(len(st), c))
    return {
        "x_nm": np.concatenate(xs),
        "copy": np.concatenate(copies),
        "start": np.concatenate(starts),
        "end": np.concatenate(ends),
    }


def _variability(copy: np.ndarray, start: np.ndarray) -> np.ndarray:
    """L'écart-type de profondeur mesuré en semaine 7, porté par locus."""
    from ..ensemble import read

    e = read(ENSEMBLE)
    sd = e.radial().std(axis=0) * e.beads.nuclear_radius * 1000.0
    out = np.empty(len(copy))
    for c in np.unique(copy):
        m = e.beads.copy_id == c
        s = e.beads.start[m]
        j = np.clip(np.searchsorted(s, start[copy == c], side="right") - 1, 0, m.sum() - 1)
        out[copy == c] = sd[m][j]
    return out


def scaling(*, cache: Path, fraction: float, sizes: str):
    """Écrit un niveau par résolution et rend les lignes du tableau."""
    bps = [int(s) for s in sizes.split(",")]
    total = genome_bp()
    yield (f"\n  {'résolution':>10} {'billes':>9} {'structure':<26} {'chunks':>6} "
           f"{'premier rendu':>14} {'niveau entier':>14} {'o/bille':>8}")
    base100 = None
    for bp in bps:
        if bp >= REAL_DOWN_TO:
            s = _nucleus(bp, cache)
            kind = "noyau S6 réel"
            if bp == REAL_DOWN_TO:
                base100 = s
        else:
            if base100 is None:
                base100 = _nucleus(REAL_DOWN_TO, cache)
            _, r_nm, _ = capacity_at(total, bp)
            s = _refine(base100, REAL_DOWN_TO // bp, r_nm)
            kind = "raffinement synthétique"
        n = len(s["copy"])
        lv = Level(
            name="banc",
            bp_per_bead=bp,
            positions_nm=s["x_nm"],
            copy=s["copy"].astype(np.int64),
            start=s["start"].astype(np.int64),
            end=s["end"].astype(np.int64),
            variability_nm=_variability(s["copy"], s["start"]),
            copies=[str(c) for c in range(int(s["copy"].max()) + 1)],
            radius_nm=[1.0] * (int(s["copy"].max()) + 1),
        )
        rep = write(cache / f"banc-{bp}.g3d", [lv], fraction=fraction)
        cols = rep["columns"]["banc"]
        first = cols["copy"]["stored"] + cols["position"]["stored"]
        whole = rep["bytes"]
        label = f"{bp // 1000} kb" if bp < 1_000_000 else f"{bp // 1_000_000} Mb"
        yield (f"  {label:>10} {n:>9,} {kind:<26} {rep['chunks']['banc']:>6} "
               f"{first / 1024:>11.1f} Kio {whole / 1024:>11.1f} Kio {whole / n:>8.2f}")
