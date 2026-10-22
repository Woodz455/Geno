"""Les allèles du locus EPHA4 en 3D : le type sauvage, le DEL-B central et ses deux témoins.

Mêmes paramètres que les régions fines des semaines 8 et 12 (8 réplicats × 150 instantanés,
mêmes graines pour tous les allèles) : ce qui sépare deux exécutions est l'allèle, rien
d'autre. Le reste de la famille de variants ne passe que par la prédiction 1D (`predict`).

Deux plans. **Pré-enregistré** (`epha4-<allèle>.npz`) : l'allèle est la région moins la
délétion. **Révisé, fenêtre constante** (`epha4-<allèle>.w.npz`) : chaque allèle garde 5 Mb, la
délétion faisant entrer la séquence en aval — voir `Locus.sites`. Le type sauvage est le même
dans les deux plans (même fenêtre, mêmes sites, mêmes graines) : une seule exécution.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..fine.polymer import Field
from ..fine.run import Fine, Setup, build, save
from . import locus as L

FIRST_SEED = 15_000


def alleles(loc: L.Locus) -> dict[str, L.Variant | None]:
    d, s, k = L.central(loc)
    return {"wt": None, "delb": d, "delbs": s, "delb-keep": k}


def path(out_dir: Path, tag: str, window: bool = False) -> Path:
    """Le type sauvage du plan révisé est celui du plan pré-enregistré tant que `extend` ne lui a
    pas ajouté de réplicats."""
    w = out_dir / f"epha4-{tag}.w.npz"
    if window and (tag != "wt" or w.exists()):
        return w
    return out_dir / f"epha4-{tag}.npz"


def run(out_dir: Path, tags: list[str] | None = None, *, replicates: int = 8, snapshots: int = 150,
        workers: int = 0, window: bool = False, log=print) -> None:
    loc = L.build()
    out_dir.mkdir(parents=True, exist_ok=True)
    for tag, v in alleles(loc).items():
        if tags and tag not in tags:
            continue
        if window and v is None:
            continue                                    # le même dans les deux plans
        sites, mids = loc.sites(v, window=window)
        region = loc.allele_region(v, window=window)
        setup = Setup(region=region, replicates=replicates, snapshots=snapshots, stride=4, relax=150_000,
                      field=Field(trunc=3.0, stiffness=1.5))
        log(f"{tag:10s} {v.name if v else 'type sauvage'} : {region.n} monomères, {sites.k} sites CTCF")
        fine = build(setup, sites=sites, first_seed=FIRST_SEED, workers=workers,
                     progress=lambda d, n: log(f"  réplicat {d}/{n}"))
        out = path(out_dir, tag, window)
        save(fine, str(out))
        side = json.loads(out.with_suffix(".json").read_text())
        side.update({
            "allele": tag, "window": window,
            "variant": None if v is None else {"name": v.name, "kind": v.kind, "start": v.start, "end": v.end,
                                               "keep": [list(k) for k in v.keep]},
            "mids_hg38": [int(x) for x in mids],
            "locus_source": L.SOURCE,
        })
        out.with_suffix(".json").write_text(json.dumps(side, ensure_ascii=False))
        log(f"  {out} — {fine.elapsed / 60:.1f} min")


def extend(out_dir: Path, tags: list[str] | None = None, *, replicates: int = 48, step: int = 8,
           workers: int = 0, log=print) -> None:
    """Ajoute des réplicats au plan révisé, `step` à la fois et allèle après allèle, jusqu'à
    `replicates` chacun — chaque tour est enregistré, les allèles avancent ensemble.

    Pourquoi : un réplicat, c'est 150 000 pas de mise en place puis 30 000 d'échantillonnage ;
    ses 150 instantanés partagent donc à peu près **une** conformation à grande échelle. À
    8 réplicats, un contact enhancer → promoteur à 1 Mb tient à un ou deux d'entre eux (erreur
    type de 50 à 100 % de la moyenne). Les graines continuent la suite : `FIRST_SEED + k`."""
    loc = L.build()
    todo = {t: v for t, v in alleles(loc).items() if not tags or t in tags}
    while True:
        progressed = False
        for tag, v in todo.items():
            src = path(out_dir, tag, window=True)
            z = np.load(str(src))
            side = json.loads(src.with_suffix(".json").read_text())
            have = len(z["seeds"])
            if have >= replicates:
                continue
            sites, mids = loc.sites(v, window=True)
            if not np.array_equal(sites.pos, z["site_pos"]):
                raise ValueError(f"{tag} : les sites ont changé depuis l'exécution enregistrée")
            region = loc.allele_region(v, window=True)
            add = min(step, replicates - have)
            setup = Setup(region=region, replicates=add, snapshots=side["snapshots"], stride=side["stride"],
                          relax=side["relax"], field=Field(trunc=3.0, stiffness=1.5))
            log(f"{tag:10s} réplicats {have + 1}–{have + add} (graines {FIRST_SEED + have}–{FIRST_SEED + have + add - 1})")
            fine = build(setup, sites=sites, first_seed=FIRST_SEED + have, workers=workers)
            total = have + add
            merged = Fine(
                coords=np.concatenate([z["coords"], fine.coords.astype(np.float32)]),
                setup=Setup(region=region, replicates=total, snapshots=side["snapshots"], stride=side["stride"],
                            relax=side["relax"], field=Field(trunc=3.0, stiffness=1.5)),
                sites=sites, truth=None, seeds=np.r_[z["seeds"], fine.seeds],
                elapsed=side["elapsed_s"] + fine.elapsed,
                lef_occupancy=(side["lef_occupancy"] * have + fine.lef_occupancy * add) / total,
                legs=np.concatenate([z["legs"].astype(np.int64), fine.legs]),
            )
            out = out_dir / f"epha4-{tag}.w.npz"
            tmp = out.with_name(out.stem + ".tmp.npz")
            save(merged, str(tmp))
            new = json.loads(tmp.with_suffix(".json").read_text())
            new.update({k: side[k] for k in ("allele", "variant", "mids_hg38", "locus_source")})
            new["window"] = True
            tmp.with_suffix(".json").write_text(json.dumps(new, ensure_ascii=False))
            tmp.replace(out)
            tmp.with_suffix(".json").replace(out.with_suffix(".json"))
            log(f"  {out.name} : {total} réplicats — {fine.elapsed / 60:.1f} min")
            progressed = True
        if not progressed:
            return


def load(file: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Coordonnées (frames, n, 3) en sigma, jambes (frames, n_lef, 2), métadonnées."""
    z = np.load(str(file))
    meta = json.loads(Path(str(file).replace(".npz", ".json")).read_text())
    return z["coords"], z["legs"], meta
