"""Produire N repliements du même génome, en parallèle, sans perdre le travail fait.

Une structure coûte une trentaine de secondes ; deux cents en coûtent une
demi-heure sur quatre cœurs. Deux conséquences sur la forme du code :

- **parallèle** — les structures sont indépendantes par construction, chacune
  part de sa propre graine et ne regarde aucune autre. C'est le cas le plus
  simple qui soit, et il serait absurde de ne pas en profiter ;
- **incrémental** — chaque structure est rangée dès qu'elle sort. Une
  interruption coûte une structure, pas la série, et `--resume` reprend là où on
  s'était arrêté.

Le génome, lui, est fixé une fois pour toutes par `lad_seed`. Chaque ouvrier
renvoie une empreinte du sien et le parent la compare : deux cents structures
d'un génome qui dérive seraient un ensemble de rien du tout, et rien dans les
coordonnées ne permettrait de s'en apercevoir après coup.
"""

from __future__ import annotations

import multiprocessing as mp
import time
from pathlib import Path
from typing import Callable

import numpy as np

from ..nucleus.beads import segment
from ..nucleus.build import build
from ..nucleus.karyotype import gm12878
from . import store

_SETTINGS: dict = {}


def _init(settings: dict) -> None:
    _SETTINGS.clear()
    _SETTINGS.update(settings)


def _one(task: tuple[int, int]):
    index, seed = task
    t0 = time.perf_counter()
    nucleus = build(seed=seed, **_SETTINGS)
    b = nucleus.beads
    return (
        index,
        seed,
        nucleus.x.astype(np.float32),
        nucleus.final,
        (b.n, round(float(b.lad.sum()), 6)),      # empreinte du génome
        time.perf_counter() - t0,
    )


def generate(
    out: Path,
    *,
    n_structures: int = 200,
    workers: int = 0,
    lad_seed: int = 0,
    first_seed: int = 1_000,
    fresh: bool = False,
    progress: Callable[[int, int, float, object], None] | None = None,
    **build_kwargs,
) -> Path:
    """Remplit `out` avec `n_structures` repliements du même génome.

    Reprend un magasin existant par défaut ; `fresh=True` l'efface et recommence.
    """
    out = Path(out)
    workers = workers or mp.cpu_count()
    karyotype = gm12878()

    settings = dict(build_kwargs)
    settings["lad_seed"] = lad_seed
    beads = segment(
        karyotype,
        bp_per_bead=settings.get("bp_per_bead", 750_000),
        nuclear_radius=settings.get("nuclear_radius", 5.0),
        phi=settings.get("phi", 0.30),
        seed=lad_seed,
    )
    fingerprint = (beads.n, round(float(beads.lad.sum()), 6))

    # Reprendre est le **défaut**, et recréer demande de le dire. Une série de deux
    # cents structures coûte une demi-heure : un `make ensemble` lancé deux fois ne
    # doit pas effacer la première, et `zarr.open_group(mode="w")` le ferait sans
    # rien demander. Un magasin déjà complet ne produit alors plus rien du tout,
    # et l'appelant enchaîne directement sur le rapport.
    if not fresh and (out / ".zgroup").exists():
        import zarr

        root = zarr.open_group(str(out), mode="a")
        stored = int(root.attrs.get("n_structures", 0))
        if stored != n_structures:
            raise ValueError(
                f"{out} contient {stored} structures, {n_structures} demandées. "
                "Passer --fresh pour le recréer, ou reprendre le même nombre."
            )
        todo = store.pending(out)
    else:
        root = store.create(
            out,
            beads,
            n_structures,
            {
                "karyotype": karyotype.name,
                "assembly": karyotype.assembly,
                "chrom_sizes_provenance": karyotype.provenance,
                "lad_seed": int(lad_seed),
                "first_seed": int(first_seed),
                "workers": int(workers),
                "build_settings": {
                    k: v for k, v in sorted(settings.items()) if not callable(v)
                },
            },
        )
        todo = np.arange(n_structures)

    tasks = [(int(i), int(first_seed + i)) for i in todo]
    if not tasks:
        return out

    started = time.perf_counter()
    with mp.Pool(workers, initializer=_init, initargs=(settings,)) as pool:
        for k, result in enumerate(pool.imap_unordered(_one, tasks), start=1):
            index, seed, coords, quality, seen, elapsed = result
            if seen != fingerprint:
                raise RuntimeError(
                    f"la structure {index} (graine {seed}) porte un autre génome : "
                    f"{seen} au lieu de {fingerprint}. L'ensemble serait ininterprétable."
                )
            store.put(root, index, seed, coords, quality)
            if progress is not None:
                progress(k, len(tasks), time.perf_counter() - started, quality)
    return out
