"""C5 et C6 : la position radiale du modèle contre la lamina (pA-DamID LMNB1, **K562, pas
GM12878**) et contre la FISH publiée des territoires en lymphoblastes.

`ensemble.stats.damid`, écrit en semaine 7, attendait une piste DamID réelle sans jamais en avoir
reçu. Il corrèle en Pearson ; les critères de la semaine 21 sont en Spearman. La lecture de la piste
est donc refaite ici, avec le lecteur bigWig de la fiche (`fiche.sources.BigWig`).
"""

from __future__ import annotations

import numpy as np

from ..fiche.sources import BigWig, locked, manifest
from .common import Model, Result, null_p95, spearman, verdict
from .sprite import CHROMS


def regions(model: Model) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Les régions haploïdes (billes de la première copie) : chromosome, et les billes des deux copies."""
    chrom, first, second = [], [], []
    for ch in CHROMS:
        cs = model.copies(ch)
        a = np.flatnonzero(model.copy_id == cs[0])
        b = np.flatnonzero(model.copy_id == cs[1]) if len(cs) > 1 else a
        n = min(len(a), len(b))
        chrom += [ch] * n; first.append(a[:n]); second.append(b[:n])
    return np.array(chrom), np.concatenate(first), np.concatenate(second)


def criteria(model: Model) -> tuple[list[Result], dict]:
    res = []
    rb = model.radial()
    chrom, a, b = regions(model)
    r_model = (rb[a] + rb[b]) / 2
    bw = BigWig(locked(manifest(), "lads.damid.lmnb1").path)
    signal = np.array([bw.summary(c, int(model.start[i]), int(model.end[i]))[0] for c, i in zip(chrom, a)])
    groups = [np.flatnonzero(chrom == ch) for ch in CHROMS]

    def permute(rng):
        p = np.arange(len(chrom))
        for s in groups:
            p[s] = rng.permutation(s)
        return r_model[p]

    r = spearman(r_model, signal)
    p95, mean = null_p95(lambda v: spearman(v, signal), permute, sign=-1)
    # Témoin : la même mesure contre la piste LAD synthétique que le noyau a reçue — une
    # corrélation imposée, que la mesure doit voir.
    lad = (model.lad[a] + model.lad[b]) / 2
    res.append(Result("C5", "pA-DamID LMNB1 (K562)", "ρ position radiale / signal lamina, par bille", r,
                      "≤ −0,3", verdict(r, -0.3, sign=-1, p95=p95), null_p95=-p95,
                      extra={"beads": int(np.isfinite(signal).sum()), "null_mean": -mean,
                             "control_synthetic_lad_rho": spearman(r_model, lad),
                             "synthetic_vs_damid_rho": spearman(lad, signal)}))

    per = {ch: float(np.mean(rb[np.isin(model.copy_id, model.copies(ch))])) for ch in CHROMS}
    ok = per["chr19"] < per["chr18"]
    res.append(Result("C6", "FISH en lymphoblastes (Croft 1999, Boyle 2001 ; texte)",
                      "position radiale moyenne chr19 − chr18 (négatif : chr19 plus central)",
                      per["chr19"] - per["chr18"], "< 0", "accord" if ok else "échec",
                      extra={"chr19": per["chr19"], "chr18": per["chr18"], "per_chrom": per}))
    return res, {"radial": r_model.tolist(), "damid": signal.tolist(), "lad": lad.tolist()}
