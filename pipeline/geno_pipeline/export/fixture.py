"""Fichier témoin pour la lecture croisée Python → TypeScript.

Écrit un `.g3d` synthétique à plusieurs chunks et plusieurs niveaux, et à côté ce que le
lecteur Python en relit. Le lecteur TypeScript doit retrouver exactement les mêmes
entiers et les mêmes positions au flottant près : deux implémentations du même format
qui ne s'accorderaient pas sur un fichier commun feraient deux formats.

    python -m geno_pipeline.export.fixture <dossier>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from .build import pose
from .g3d import Level, read, write


def _chain(counts, name, bp, seed, step_nm):
    rng = np.random.default_rng(seed)
    pos, copy, start = [], [], []
    for c, k in enumerate(counts):
        pos.append(np.cumsum(rng.normal(0, step_nm, (k, 3)), axis=0) + rng.normal(0, 3000, 3))
        copy.append(np.full(k, c))
        start.append(np.arange(k) * bp + 10_000)
    n = sum(counts)
    start = np.concatenate(start)
    return Level(
        name=name,
        bp_per_bead=bp,
        positions_nm=np.concatenate(pos),
        copy=np.concatenate(copy),
        start=start,
        end=start + bp,
        variability_nm=rng.uniform(50, 500, n),
        copies=[f"chr{c // 2 + 1}:{'ab'[c % 2]}" for c in range(len(counts))],
        radius_nm=[150.0] * len(counts),
    )


def main(out: str) -> None:
    d = Path(out)
    d.mkdir(parents=True, exist_ok=True)
    nucleus = _chain((900, 900, 700, 700), "noyau", 750_000, 1, 230.0)
    finer = _chain((1500,), "fin", 2_000, 2, 30.0)
    anchors = finer.positions_nm[::300][:5]
    m, fit = pose(anchors @ np.diag([1, -1, 1]) + 500.0, anchors)
    finer.transform = m.reshape(-1).tolist()
    finer.fit = fit
    rep = write(d / "fixture.g3d", [nucleus, finer], first="noyau", max_beads=500,
                header={"assembly": "GRCh38", "warnings": ["témoin synthétique"]})
    got = read(d / "fixture.g3d")

    expected = {"bytes": rep["bytes"], "first_frame_end": rep["first_frame_end"], "levels": {}}
    for name, lv in got.levels.items():
        meta = lv["meta"]
        hom = np.c_[lv["positions_nm"], np.ones(len(lv["positions_nm"]))]
        nucleus_frame = (np.asarray(meta["frame"]["transform"]).reshape(4, 4) @ hom.T).T[:, :3]
        # Requêtes de référence, résolues par force brute sur les colonnes relues.
        queries = []
        for chrom, s, e in (("chr1", 30_000_000, 34_000_000), ("chr2", 0, 760_000),
                            ("chr1", 13_000, 14_000), ("chr9", 0, 10**9)):
            hits = []
            for c, label in enumerate(meta["copies"]):
                if label.split(":")[0] != chrom:
                    continue
                ids = np.flatnonzero((lv["copy"] == c) & (lv["start"] < e) & (lv["end"] > s))
                ids = ids[np.argsort(lv["start"][ids])]
                if len(ids):
                    hits.append({"copy": label, "ids": (ids + 1).tolist()})
            queries.append({"chrom": chrom, "start": s, "end": e, "hits": hits})
        expected["levels"][name] = {
            "n": len(lv["copy"]),
            "chunks": len(meta["chunks"]),
            "positions": lv["positions_nm"].reshape(-1).tolist(),
            "nucleus_frame": nucleus_frame.reshape(-1).tolist(),
            "copy": lv["copy"].tolist(),
            "start": lv["start"].tolist(),
            "end": lv["end"].tolist(),
            "variability": lv["variability_nm"].tolist(),
            "variability_step": meta["variability"]["step"],
            "queries": queries,
        }
    (d / "expected.json").write_text(json.dumps(expected))
    print(f"{d / 'fixture.g3d'}  {rep['bytes']:,} o  ·  chunks {rep['chunks']}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "dist/g3d-fixture")
