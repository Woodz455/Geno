"""Assemble les quatre étages et produit un ensemble de conformations fines.

Un **ensemble**, pas une structure : le principe n° 1 du projet ne change pas
d'échelle. Une carte de contacts issue d'une seule trajectoire mesurerait une
trajectoire, et la semaine 7 a chiffré ce que ça coûte au niveau du noyau.

Les réplicats sont indépendants — graine, conformation initiale, tirage
d'extrusion — et parallélisés. Ce qu'ils partagent est ce qui doit l'être : la
région, les sites CTCF, la vérité plantée. Exactement la séparation
`lad_seed`/`seed` de la semaine 7, pour la même raison : si chaque réplicat
tirait ses propres sites, la variabilité mesurée mélangerait la variabilité de
repliement et celle du génome.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import time
from dataclasses import dataclass, replace

import numpy as np

from .extrusion import simulate as extrude
from .polymer import Field
from .polymer import simulate as fold
from .region import Region, Sites, Truth, blob_radius, plant

FORMAT = "geno-fine/1"


@dataclass(frozen=True)
class Setup:
    """Tout ce qui définit une exécution, et rien d'autre."""

    region: Region
    separation: int = 200_000
    processivity: int = 200_000
    release: float = 0.003
    velocity: float = 1.0
    phi: float = 0.30
    field: Field = Field()
    md_per_step: int = 200
    relax: int = 150_000
    snapshots: int = 120
    stride: int = 4
    burnin_1d: int = 2_000
    replicates: int = 4
    cutoff: float = 1.5

    @property
    def confine(self) -> float:
        return blob_radius(self.region.n, self.phi)

    @property
    def sigma_nm(self) -> float:
        return 2.0 * self.region.radius_nm(self.phi)

    def as_dict(self) -> dict:
        return {
            "format": FORMAT,
            "chrom": self.region.chrom,
            "start": self.region.start,
            "end": self.region.end,
            "bp_per_bead": self.region.bp_per_bead,
            "monomers": self.region.n,
            "sigma_nm": round(self.sigma_nm, 3),
            "confine_sigma": round(self.confine, 3),
            "phi": self.phi,
            "separation_bp": self.separation,
            "processivity_bp": self.processivity,
            "release": self.release,
            "md_per_step": self.md_per_step,
            "relax": self.relax,
            "snapshots": self.snapshots,
            "stride": self.stride,
            "replicates": self.replicates,
            "cutoff": self.cutoff,
            "trunc_kT": self.field.trunc,
            "stiffness_kT": self.field.stiffness,
            "evidence": "simulated",
        }


@dataclass(frozen=True)
class Fine:
    """Le résultat : des conformations, et de quoi savoir d'où elles viennent."""

    coords: np.ndarray        # (frames, n, 3) en unités de sigma
    setup: Setup
    sites: Sites
    truth: Truth | None
    seeds: np.ndarray         # (replicates,)
    elapsed: float
    lef_occupancy: float

    @property
    def frames(self) -> int:
        return len(self.coords)


_SHARED: dict = {}


def _init(setup: Setup, sites: Sites) -> None:
    _SHARED["setup"] = setup
    _SHARED["sites"] = sites


def _one(args: tuple[int, int]) -> tuple[int, np.ndarray, float]:
    index, seed = args
    setup: Setup = _SHARED["setup"]
    sites: Sites = _SHARED["sites"]
    ex = extrude(
        setup.region,
        sites,
        separation=setup.separation,
        processivity=setup.processivity,
        velocity=setup.velocity,
        release=setup.release,
        burnin=setup.burnin_1d,
        snapshots=setup.snapshots,
        stride=setup.stride,
        seed=seed,
    )
    coords = fold(
        setup.region.n,
        ex.legs,
        confine=setup.confine,
        field=setup.field,
        md_per_step=setup.md_per_step,
        relax=setup.relax,
        seed=seed,
        threads=1,
    )
    return index, coords, ex.occupancy


def build(
    setup: Setup,
    *,
    sites: Sites | None = None,
    truth: Truth | None = None,
    site_seed: int = 1,
    first_seed: int = 4_000,
    workers: int = 0,
    progress=None,
    **plant_kw,
) -> Fine:
    """Fait tourner `setup.replicates` trajectoires indépendantes et les empile.

    `sites` non fourni ⇒ on plante une région de vérité connue, témoins compris.
    La graine des sites est **distincte** de celles des réplicats, et c'est le
    point : tous les réplicats voient les mêmes barrières.
    """
    if sites is None:
        sites, truth = plant(setup.region, seed=site_seed, **plant_kw)

    seeds = first_seed + np.arange(setup.replicates, dtype=np.int64)
    jobs = list(enumerate(int(s) for s in seeds))
    workers = workers or min(setup.replicates, mp.cpu_count())

    started = time.time()
    chunks: list[np.ndarray | None] = [None] * setup.replicates
    occ: list[float] = []

    if workers <= 1:
        _init(setup, sites)
        results = (_one(j) for j in jobs)
        for done, (index, coords, o) in enumerate(results, 1):
            chunks[index] = coords
            occ.append(o)
            if progress:
                progress(done, len(jobs))
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(workers, initializer=_init, initargs=(setup, sites)) as pool:
            for done, (index, coords, o) in enumerate(pool.imap_unordered(_one, jobs), 1):
                chunks[index] = coords
                occ.append(o)
                if progress:
                    progress(done, len(jobs))

    return Fine(
        coords=np.concatenate(chunks, axis=0),
        setup=setup,
        sites=sites,
        truth=truth,
        seeds=seeds,
        elapsed=time.time() - started,
        lef_occupancy=float(np.mean(occ)) if occ else 0.0,
    )


def knockout(setup: Setup, sites: Sites, **kw) -> Fine:
    """La même chose avec les sites inoccupés — le témoin négatif du modèle entier.

    Rien d'autre ne change : mêmes graines, mêmes positions de sites, même champ
    de force. Ce qui sépare les deux exécutions est donc exactement l'arrêt de
    l'extrusion par CTCF, ce qui fait de la différence entre elles une mesure et
    non une illustration.
    """
    return build(setup, sites=sites.without_ctcf(), **kw)


def save(fine: Fine, path: str) -> None:
    """`.npz` pour les tableaux, `.json` à côté pour la provenance.

    Le sidecar porte `evidence: "simulated"` comme partout ailleurs dans le
    projet, et la source des sites CTCF : `planted:<graine>` n'est pas un fichier
    publié, et rien dans la sortie ne doit laisser croire le contraire.
    """
    np.savez_compressed(
        path,
        coords=fine.coords.astype(np.float32),
        site_pos=fine.sites.pos,
        site_strand=fine.sites.strand,
        site_occupancy=fine.sites.occupancy,
        seeds=fine.seeds,
        **(
            {
                "truth_left": fine.truth.left,
                "truth_right": fine.truth.right,
                "truth_kind": fine.truth.kind,
                "truth_boundary": fine.truth.boundary,
            }
            if fine.truth is not None
            else {}
        ),
    )
    meta = fine.setup.as_dict()
    meta.update(
        {
            "frames": fine.frames,
            "seeds": [int(s) for s in fine.seeds],
            "ctcf_source": fine.sites.source,
            "lef_occupancy": round(fine.lef_occupancy, 4),
            "elapsed_s": round(fine.elapsed, 1),
        }
    )
    with open(str(path).replace(".npz", "") + ".json", "w") as fh:
        json.dump(meta, fh, indent=2)


def load(path: str) -> tuple[np.ndarray, Sites, Truth | None, dict]:
    z = np.load(path)
    with open(str(path).replace(".npz", "") + ".json") as fh:
        meta = json.load(fh)
    sites = Sites(
        z["site_pos"], z["site_strand"], z["site_occupancy"], meta.get("ctcf_source", "?")
    )
    truth = (
        Truth(
            z["truth_left"],
            z["truth_right"],
            z["truth_kind"],
            z["truth_boundary"],
        )
        if "truth_left" in z
        else None
    )
    return z["coords"], sites, truth, meta


def region_of(meta: dict) -> Region:
    return Region(meta["chrom"], meta["start"], meta["end"], meta["bp_per_bead"])


def scaled(setup: Setup, **changes) -> Setup:
    return replace(setup, **changes)
