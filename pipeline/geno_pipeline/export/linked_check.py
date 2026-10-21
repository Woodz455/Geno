"""Ce que vaut la carte de contacts des vues liées, mesuré par un autre chemin.

1. **Conservation** — à chaque résolution, la somme des tuiles (miroir compris) égale la
   somme de la matrice de bins : la répartition sur la grille ne perd ni ne crée de contact.
2. **Exactitude** — des pixels tirés au hasard, recomposés **sans KD-tree ni matrice
   creuse** : les distances des paires de billes qu'ils couvrent, calculées directement sur
   les 200 repliements, comptées, réparties au prorata des recouvrements. Témoin : le même
   calcul décalé d'un pixel doit diverger.
"""

from __future__ import annotations

import random

import numpy as np

from .linked import TILE, Contacts, genome

# Les tuiles sont en float16 (11 bits de mantisse : 4,9·10⁻⁴ d'erreur relative au pire),
# plus un plancher absolu pour les sous-normaux.
TOL = 1e-3


def conservation(c: Contacts) -> list[dict]:
    out = []
    for ts in c.header["tilesets"]:
        for lv in ts["levels"]:
            s = 0.0
            for key in lv["tiles"]:
                z, x, y = map(int, key.split("."))
                t = c.tile(ts["uid"], z, x, y)
                v = np.nan_to_num(t, nan=0.0).astype(np.float64)
                if x == y:
                    s += v.sum()
                else:
                    s += 2.0 * v.sum()           # la tuile et son miroir
            out.append({"tileset": ts["uid"], "resolution": lv["resolution"], "tiles": len(lv["tiles"]),
                        "stored": s, "expected": lv["total"], "rel": abs(s - lv["total"]) / max(lv["total"], 1e-12)})
    return out


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def exactness(c: Contacts, ensemble_path, n: int = 200, seed: int = 14, shift: int = 0,
              cutoff: float = 1.5) -> dict:
    """Pixels de 750 kb : un tiers sur la diagonale, un tiers proches (cis), un tiers trans."""
    from ..ensemble.compartments import haploid_bins
    from ..ensemble.store import read

    e = read(ensemble_path)
    b = e.beads
    sizes = dict(c.header["chromsizes"])
    offsets, total = genome(sizes)
    bin_of, bins, _ = haploid_bins(tuple(b.labels), b.copy_id, b.start, b.end)
    off = np.array([offsets[ch] for ch in bins["chrom"]], dtype=np.int64)
    a_start = off + bins["start"].to_numpy(np.int64)
    a_end = off + bins["end"].to_numpy(np.int64)
    beads_of = [np.flatnonzero(bin_of == k) for k in range(len(bins))]
    ts = c.tileset("noyau")
    z = ts["resolutions"].index(750_000)
    res = 750_000
    n_px = -(-total // res)
    rng = random.Random(seed)
    picks = []
    while len(picks) < n:
        kind = len(picks) % 3
        p = rng.randrange(n_px)
        q = p if kind == 0 else min(n_px - 1, p + rng.randrange(1, 8)) if kind == 1 else rng.randrange(n_px)
        picks.append((max(p, q), min(p, q)))
    good = 0
    first = None
    for row, col in picks:
        def bins_in(px):
            lo, hi = px * res, (px + 1) * res
            ks = np.flatnonzero((a_start < hi) & (a_end > lo))
            return [(k, _overlap(a_start[k], a_end[k], lo, hi) / (a_end[k] - a_start[k])) for k in ks]

        want = 0.0
        for ka, fa in bins_in(row):
            for kb, fb in bins_in(col):
                # Paires de billes (diploïdes) des deux bins, comptées une fois par paire non
                # ordonnée, comme `count_contacts` (bin1 ≤ bin2).
                ia, ib = beads_of[ka], beads_of[kb]
                d = np.linalg.norm(e.coords[:, ia, None, :] - e.coords[:, None, ib, :], axis=3)
                lim = cutoff * (b.radius[ia, None] + b.radius[None, ib])
                hit = d < lim
                if ka == kb:
                    hit = hit & np.triu(np.ones(hit.shape[1:], bool), k=1)
                want += hit.sum() / e.n_structures * fa * fb
        got = c.pixel("noyau", z, row + shift, col + shift) if row + shift < n_px else 0.0
        ok = abs(got - want) <= TOL * abs(want) + 1e-7
        good += ok
        if not ok and first is None:
            first = {"row": row, "col": col, "got": got, "want": want}
    return {"pixels": n, "equal": good, "first_diff": first, "shift": shift}


def fine_exactness(c: Contacts, npz_paths, n: int = 120, seed: int = 15, shift: int = 0) -> dict:
    """Pixels de 2 kb des régions fines : une paire de monomères = un pixel (à la répartition
    près, la région ne commençant pas sur un multiple de 2 kb absolu) — distances directes."""
    from ..fine.run import load

    ts = c.tileset("fin")
    z = ts["resolutions"].index(2_000)
    res = 2_000
    rng = random.Random(seed)
    good = total = 0
    first = None
    for path in npz_paths:
        coords, _, _, meta = load(str(path))
        reg = next(r for r in ts["regions"] if r["chrom"] == meta["chrom"] and r["start"] == meta["start"])
        cutoff = float(meta.get("cutoff", 1.5))
        bp = int(meta["bp_per_bead"])
        m = coords.shape[1]
        for _ in range(n // max(1, len(npz_paths))):
            px_r = (reg["abs_start"] + rng.randrange(m * bp)) // res
            px_c = max(reg["abs_start"] // res, px_r - rng.randrange(0, 60))
            want = 0.0

            def monomers(px):
                lo, hi = px * res, (px + 1) * res
                first = max(0, (lo - reg["abs_start"]) // bp)
                last = min(m - 1, (hi - 1 - reg["abs_start"]) // bp)
                for i in range(first, last + 1):
                    s0 = reg["abs_start"] + i * bp
                    f = _overlap(s0, s0 + bp, lo, hi) / bp
                    if f:
                        yield i, f

            for i, fi in monomers(px_r):
                for j, fj in monomers(px_c):
                    if i != j:                           # un monomère n'est pas en contact avec lui-même
                        d = np.linalg.norm(coords[:, i] - coords[:, j], axis=1)
                        want += (d < cutoff).mean() * fi * fj
            got = c.pixel("fin", z, px_r + shift, px_c + shift)
            got = 0.0 if np.isnan(got) else got
            ok = abs(got - want) <= TOL * abs(want) + 1e-7
            good += ok
            total += 1
            if not ok and first is None:
                first = {"region": reg["name"], "row": px_r, "col": px_c, "got": got, "want": want}
    return {"pixels": total, "equal": good, "first_diff": first, "shift": shift}


def dump_expected(c: Contacts, out, n: int = 40, seed: int = 16) -> None:
    """Ce que lit le lecteur Python, pour la lecture croisée TypeScript (`tools/test-contacts.mjs`) :
    l'empreinte du float32 décodé de tuiles tirées au hasard, et des pixels — miroir et tuiles
    absentes compris."""
    import hashlib
    import json
    from pathlib import Path

    rng = random.Random(seed)
    tiles, pixels = [], []
    for ts in c.header["tilesets"]:
        keys = [k for lv in ts["levels"] for k in lv["tiles"]]
        for key in sorted(rng.sample(keys, min(n, len(keys)))):
            z, x, y = map(int, key.split("."))
            t = c.tile(ts["uid"], z, x, y)
            tiles.append({"uid": ts["uid"], "z": z, "x": x, "y": y,
                          "sha256": hashlib.sha256(np.ascontiguousarray(t, "<f4").tobytes()).hexdigest(),
                          "nan": int(np.isnan(t).sum())})
        for _ in range(n):
            lv = rng.choice(ts["levels"])
            r, q = rng.randrange(lv["pixels"]), rng.randrange(lv["pixels"])
            v = c.pixel(ts["uid"], lv["z"], r, q)
            pixels.append({"uid": ts["uid"], "z": lv["z"], "row": r, "col": q, "value": None if np.isnan(v) else v})
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps({"tiles": tiles, "pixels": pixels}))


__all__ = ["conservation", "exactness", "fine_exactness", "dump_expected", "TILE"]
