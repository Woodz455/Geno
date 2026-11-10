"""Un `.mcool` à structure connue : ce qu'on importe quand on veut pouvoir juger le résultat.

Deux chromosomes fictifs, chacun une conformation de la semaine 5 (marche persistante confinée,
compartiments, volume exclu), ses contacts tirés par le modèle direct `f ∝ d^(−gamma)` avec un
biais de couverture par bin — puis un vrai fichier cooler multi-résolution, équilibré comme le
sont ceux de 4DN. La vérité (les coordonnées) est écrite à côté ; la session ne la lit jamais.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..hic.polymer import chain, contacts

CHROMS = ("chrP1", "chrP2")


def plant(out: Path, *, n: tuple[int, int] = (600, 400), resolution: int = 50_000, gamma: float = 3.0,
          total: int = 40_000_000, bias_sd: float = 0.3, seed: int = 0, factors: tuple[int, ...] = (1, 2, 4)) -> dict:
    import cooler

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    bins, pix, truth, offset = [], [], {}, 0
    for k, (name, nk) in enumerate(zip(CHROMS, n)):
        conf = chain(n=nk, seed=seed + k)
        c = contacts(conf, gamma=gamma, total=total * nk // sum(n), seed=seed + 10 + k, bias_sd=bias_sd)
        i, j = np.nonzero(np.triu(c, 1))
        pix.append(pd.DataFrame({"bin1_id": i + offset, "bin2_id": j + offset, "count": c[i, j].astype(np.int32)}))
        bins.append(pd.DataFrame({"chrom": name, "start": np.arange(nk) * resolution, "end": (np.arange(nk) + 1) * resolution}))
        truth[name] = conf.coords
        offset += nk
    bins = pd.concat(bins, ignore_index=True)
    pixels = pd.concat(pix, ignore_index=True).sort_values(["bin1_id", "bin2_id"]).reset_index(drop=True)
    base = out.with_suffix(".base.cool")
    cooler.create_cooler(str(base), bins=bins, pixels=pixels, dtypes={"count": "int32"}, assembly="planted",
                         ordered=True, symmetric_upper=True)
    if out.exists():
        out.unlink()
    res = [resolution * f for f in factors]
    cooler.zoomify_cooler(str(base), str(out), res, chunksize=10_000_000)
    base.unlink()
    for r in res:
        clr = cooler.Cooler(f"{out}::resolutions/{r}")
        cooler.balance_cooler(clr, store=True, ignore_diags=1, min_nnz=5, mad_max=0, cis_only=True)
    side = out.with_suffix(".truth.npz")
    np.savez(side, **truth)
    meta = {"chroms": {c: int(k) for c, k in zip(CHROMS, n)}, "resolution": resolution, "resolutions": res,
            "gamma": gamma, "total": total, "bias_sd": bias_sd, "seed": seed}
    out.with_suffix(".truth.json").write_text(json.dumps(meta))
    return meta


def truth(mcool: Path, chrom: str, factor: int = 1) -> np.ndarray:
    """La vérité à la résolution de base × `factor` : moyenne des billes de chaque bin."""
    x = np.load(Path(mcool).with_suffix(".truth.npz"))[chrom]
    n = len(x) // factor * factor
    out = x[:n].reshape(-1, factor, 3).mean(axis=1)
    return np.vstack([out, x[n:].mean(axis=0, keepdims=True)]) if n < len(x) else out
