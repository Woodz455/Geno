"""`ensemble.zarr` — N structures, le génome qu'elles partagent, et d'où tout vient.

Un ensemble n'est pas une pile de fichiers de structures. C'est **un génome** et
**N repliements** de ce génome : le découpage en billes, les rayons et la piste
LAD sont écrits **une seule fois**, et seules les coordonnées portent l'indice de
structure. Ce n'est pas une économie de place — c'est ce qui rend l'objet
interprétable. Si chaque structure portait son propre génome, la variabilité
mesurée mélangerait deux choses, et il n'y aurait aucun moyen de les séparer
après coup.

Écriture **incrémentale** : les tableaux sont alloués d'emblée et remplis au fur
et à mesure, avec un masque `done`. Une génération de 200 structures dure une
demi-heure ; elle doit survivre à son interruption, et pouvoir reprendre.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import zarr

from ..nucleus.beads import Beads

FORMAT = "geno-ensemble/1"


@dataclass(frozen=True)
class Ensemble:
    """Un ensemble chargé : un génome, N repliements, et leur qualité."""

    path: Path
    beads: Beads
    coords: np.ndarray        # (N, n, 3) µm
    seeds: np.ndarray         # (N,)
    max_overlap: np.ndarray   # (N,)
    bond_stretch: np.ndarray  # (N,)
    outside: np.ndarray       # (N,)
    shakes: np.ndarray        # (N,)
    meta: dict

    @property
    def n_structures(self) -> int:
        return len(self.coords)

    @property
    def n_beads(self) -> int:
        return self.beads.n

    def radial(self) -> np.ndarray:
        """(N, n) position radiale normalisée — **la** quantité comparable entre
        structures, parce qu'elle ne dépend pas de l'orientation (§ `stats`)."""
        return (
            np.linalg.norm(self.coords, axis=2) / self.beads.nuclear_radius
        ).astype(np.float64)

    def __str__(self) -> str:
        return (
            f"{self.n_structures} structures × {self.n_beads:,} billes · "
            f"{self.meta.get('karyotype', '?')} · "
            f"piste LAD {self.meta.get('lad_source', '?')} · "
            f"longueurs {self.meta.get('chrom_sizes_provenance', '?')}"
        )


def create(path: Path, beads: Beads, n_structures: int, meta: dict) -> zarr.Group:
    """Alloue le magasin et y écrit le génome partagé. Les coordonnées viennent après."""
    path = Path(path)
    root = zarr.open_group(str(path), mode="w")

    root.attrs.update(
        {
            "format": FORMAT,
            "evidence": "simulated",
            "n_structures": int(n_structures),
            "n_beads": int(beads.n),
            "bp_per_bead": int(np.median(beads.end - beads.start)),
            "nuclear_radius_um": float(beads.nuclear_radius),
            "phi": float(beads.phi),
            "mean_bead_radius_nm": round(float(beads.radius.mean()) * 1000, 1),
            "lad_source": beads.lad_source,
            "labels": list(beads.labels),
            "warning": (
                "Ensemble simulé. Aucune donnée de conformation ne contraint ces "
                "positions ; la piste LAD est synthétique. Une structure isolée de "
                "cet ensemble n'est pas « le » génome — c'est un tirage."
            ),
            **meta,
        }
    )

    g = root.create_group("beads")
    for name, arr in (
        ("copy_id", beads.copy_id),
        ("start", beads.start),
        ("end", beads.end),
        ("radius", beads.radius.astype(np.float32)),
        ("lad", beads.lad.astype(np.float32)),
        ("acrocentric", beads.acrocentric),
    ):
        g.create_dataset(name, data=np.asarray(arr), overwrite=True)

    # Une structure par chunk : on écrit et on relit structure par structure.
    root.create_dataset(
        "coords",
        shape=(n_structures, beads.n, 3),
        chunks=(1, beads.n, 3),
        dtype="f4",
        overwrite=True,
    )
    root.create_dataset("done", data=np.zeros(n_structures, dtype=bool), overwrite=True)
    q = root.create_group("quality")
    for name, dtype in (
        ("seed", "i8"),
        ("max_overlap", "f4"),
        ("bond_stretch", "f4"),
        ("outside", "i4"),
        ("shakes", "i4"),
    ):
        q.create_dataset(name, shape=(n_structures,), dtype=dtype, overwrite=True)
    return root


def put(root: zarr.Group, index: int, seed: int, coords: np.ndarray, quality) -> None:
    """Range une structure. `done` est écrit **en dernier** : un magasin interrompu
    ne prétend jamais contenir une structure à moitié écrite."""
    root["coords"][index] = coords.astype(np.float32)
    q = root["quality"]
    q["seed"][index] = seed
    q["max_overlap"][index] = quality.max_overlap
    q["bond_stretch"][index] = quality.bond_stretch
    q["outside"][index] = quality.outside
    q["shakes"][index] = quality.shakes
    root["done"][index] = True


def read(path: Path) -> Ensemble:
    """Relit un ensemble, **sans les structures inachevées**."""
    path = Path(path)
    root = zarr.open_group(str(path), mode="r")
    if root.attrs.get("format") != FORMAT:
        raise ValueError(f"{path} n'est pas un ensemble Geno ({root.attrs.get('format')!r})")

    done = np.asarray(root["done"][:], dtype=bool)
    if not done.any():
        raise ValueError(f"{path} ne contient aucune structure terminée")

    g = root["beads"]
    beads = Beads(
        copy_id=np.asarray(g["copy_id"][:]),
        start=np.asarray(g["start"][:]),
        end=np.asarray(g["end"][:]),
        radius=np.asarray(g["radius"][:], dtype=np.float64),
        lad=np.asarray(g["lad"][:], dtype=np.float64),
        labels=tuple(root.attrs["labels"]),
        acrocentric=np.asarray(g["acrocentric"][:], dtype=bool),
        nuclear_radius=float(root.attrs["nuclear_radius_um"]),
        phi=float(root.attrs["phi"]),
        lad_source=str(root.attrs["lad_source"]),
    )

    q = root["quality"]
    keep = np.flatnonzero(done)
    return Ensemble(
        path=path,
        beads=beads,
        coords=np.asarray(root["coords"][:])[keep].astype(np.float64),
        seeds=np.asarray(q["seed"][:])[keep],
        max_overlap=np.asarray(q["max_overlap"][:])[keep],
        bond_stretch=np.asarray(q["bond_stretch"][:])[keep],
        outside=np.asarray(q["outside"][:])[keep],
        shakes=np.asarray(q["shakes"][:])[keep],
        meta=dict(root.attrs),
    )


def pending(path: Path) -> np.ndarray:
    """Indices restant à produire dans un magasin existant, pour reprendre."""
    root = zarr.open_group(str(Path(path)), mode="r")
    return np.flatnonzero(~np.asarray(root["done"][:], dtype=bool))
