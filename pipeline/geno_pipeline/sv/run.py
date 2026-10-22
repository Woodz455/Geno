"""Les allèles du locus EPHA4 en 3D : le type sauvage, le DEL-B central et ses deux témoins.

Mêmes paramètres que les régions fines des semaines 8 et 12 (8 réplicats × 150 instantanés,
mêmes graines pour tous les allèles) : ce qui sépare deux exécutions est l'allèle, rien
d'autre. Le reste de la famille de variants ne passe que par la prédiction 1D (`predict`).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..fine.polymer import Field
from ..fine.run import Setup, build, save
from . import locus as L

FIRST_SEED = 15_000


def alleles(loc: L.Locus) -> dict[str, L.Variant | None]:
    d, s, k = L.central(loc)
    return {"wt": None, "delb": d, "delbs": s, "delb-keep": k}


def run(out_dir: Path, tags: list[str] | None = None, *, replicates: int = 8, snapshots: int = 150,
        workers: int = 0, log=print) -> None:
    loc = L.build()
    out_dir.mkdir(parents=True, exist_ok=True)
    for tag, v in alleles(loc).items():
        if tags and tag not in tags:
            continue
        sites, mids = loc.sites(v)
        region = loc.allele_region(v)
        setup = Setup(region=region, replicates=replicates, snapshots=snapshots, stride=4, relax=150_000,
                      field=Field(trunc=3.0, stiffness=1.5))
        log(f"{tag:10s} {v.name if v else 'type sauvage'} : {region.n} monomères, {sites.k} sites CTCF")
        fine = build(setup, sites=sites, first_seed=FIRST_SEED, workers=workers,
                     progress=lambda d, n: log(f"  réplicat {d}/{n}"))
        path = out_dir / f"epha4-{tag}.npz"
        save(fine, str(path))
        side = json.loads(path.with_suffix(".json").read_text())
        side.update({
            "allele": tag,
            "variant": None if v is None else {"name": v.name, "kind": v.kind, "start": v.start, "end": v.end,
                                               "keep": [list(k) for k in v.keep]},
            "mids_hg38": [int(x) for x in mids],
            "locus_source": L.SOURCE,
        })
        path.with_suffix(".json").write_text(json.dumps(side, ensure_ascii=False))
        log(f"  {path} — {fine.elapsed / 60:.1f} min")


def load(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Coordonnées (frames, n, 3) en sigma, jambes (frames, n_lef, 2), métadonnées."""
    z = np.load(str(path))
    meta = json.loads(Path(str(path).replace(".npz", ".json")).read_text())
    return z["coords"], z["legs"], meta
