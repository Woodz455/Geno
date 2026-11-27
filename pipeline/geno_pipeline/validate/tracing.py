"""C3 et C4 : le modèle contre le traçage multiplexé de Su et al. 2020 — **IMR90, pas GM12878**.

- `chromosome21.tsv` : 651 segments de 50 kb de chr21, copie par copie (7 591 copies) ;
- `genomic-scale.tsv` : 1 041 loci sur tout le génome, par homologue et par cellule (1 787
  cellules, trois expériences).

Les coordonnées sont en hg38, en nm. Une valeur manquante (spot non détecté) est NaN, jamais 0.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..ensemble.stats import slope
from ..fiche.sources import locked, manifest
from .common import Model, Result, null_p95, pair_distances, spearman, verdict
from .sprite import CHROMS, oe_vec


def _mid(coord: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    parts = coord.str.extract(r"^(chr[0-9XY]+):(\d+)-(\d+)$")
    s, e = parts[1].astype(np.int64), parts[2].astype(np.int64)
    return parts[0].to_numpy(), ((s - 1 + e) // 2).to_numpy()


def _block_median(T: np.ndarray, rows: np.ndarray, block: int = 6) -> np.ndarray:
    """Distance médiane (µm) entre chaque locus de `rows` et tous les loci, sur l'axe 0 de `T`
    (copies, ou cellules × homologues) ; `T` : (échantillons, loci, [h,] 3) en nm."""
    n = T.shape[1]
    out = np.full((len(rows), n), np.nan, dtype=np.float32)
    for b0 in range(0, len(rows), block):
        r = rows[b0:b0 + block]
        if T.ndim == 3:
            d = np.linalg.norm(T[:, r, None, :] - T[:, None, :, :], axis=3)
        else:   # (cellules, loci, homologues, 3) : toutes les paires d'homologues
            d = np.concatenate([np.linalg.norm(T[:, r, ha, None, :] - T[:, None, :, hb, :], axis=3)
                                for ha in range(T.shape[2]) for hb in range(T.shape[2])], axis=0)
        out[b0:b0 + len(r)] = np.nanmedian(d, axis=0) / 1000.0
    return out


def _nu(sep: np.ndarray, dist: np.ndarray, lo: float = 7.5e5, hi: float = 1e7) -> tuple[float, list, list]:
    """Exposant de R(s) ∝ s^ν entre `lo` et `hi`, sur des bins logarithmiques de séparation."""
    ok = np.isfinite(dist) & (sep > 0)
    edges = np.logspace(np.log10(lo), np.log10(hi), 13)
    xs, ys = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = ok & (sep >= a) & (sep < b)
        if sel.sum() >= 3:
            xs.append(np.sqrt(a * b)); ys.append(float(np.median(dist[sel])))
    nu = slope(np.array(xs), np.array(ys), lo, hi) if len(xs) >= 3 else float("nan")
    return nu, xs, ys


@dataclass
class Chr21:
    pos: np.ndarray        # (651,) centres hg38
    med: np.ndarray        # (651, 651) µm, traçage
    copies: int


def load_chr21() -> Chr21:
    df = pd.read_csv(locked(manifest(), "su2020.chr21").path, sep="\t",
                     usecols=["Z(nm)", "X(nm)", "Y(nm)", "Genomic coordinate", "Chromosome copy number"])
    _, mid = _mid(df["Genomic coordinate"])
    pos = np.unique(mid)
    li = np.searchsorted(pos, mid)
    cp, ci = np.unique(df["Chromosome copy number"].to_numpy(), return_inverse=True)
    T = np.full((len(cp), len(pos), 3), np.nan, dtype=np.float32)
    T[ci, li] = df[["Z(nm)", "X(nm)", "Y(nm)"]].to_numpy(dtype=np.float32)
    med = _block_median(T, np.arange(len(pos)))
    return Chr21(pos, med, len(cp))


def chr21_criteria(model: Model, t: Chr21) -> list[Result]:
    res = []
    idx = model.beads_at("chr21", t.pos)             # (2, 651)
    ca, cb = model.copies("chr21")
    ids_a, ids_b = np.flatnonzero(model.copy_id == ca), np.flatnonzero(model.copy_id == cb)
    beads = np.unique(idx[0][idx[0] >= 0])            # billes de la copie a qui portent un locus
    la = np.searchsorted(ids_a, beads)
    bidx = np.array([ids_a[la], ids_b[la]])           # les mêmes billes sur la copie b
    D = pair_distances(model, bidx, cis_only=True)    # (2S, B, B)
    mb = np.nanmedian(D, axis=0)
    mid_b = (model.start[beads] + model.end[beads]) / 2

    # Le modèle au niveau des loci : la distance de leurs billes, NaN dans la même bille.
    li = np.searchsorted(beads, idx[0])
    valid = idx[0] >= 0
    ml = np.where(valid[:, None] & valid[None, :], mb[np.clip(li, 0, len(beads) - 1)][:, np.clip(li, 0, len(beads) - 1)], np.nan)
    ml[li[:, None] == li[None, :]] = np.nan

    iu, ju = np.triu_indices(len(t.pos), 1)
    sep = (t.pos[ju] - t.pos[iu]).astype(float)
    nu_t, xt, yt = _nu(sep, t.med[iu, ju])
    bi, bj = np.triu_indices(len(beads), 1)
    nu_m, xm, ym = _nu((mid_b[bj] - mid_b[bi]).astype(float), mb[bi, bj])
    res.append(Result("C3a", "Su 2020 chr21 (IMR90)", "exposant ν de R(s), 0,75–10 Mb : modèle − traçage",
                      nu_m - nu_t, "|écart| ≤ 0,1", "accord" if abs(nu_m - nu_t) <= 0.1 else "échec",
                      extra={"model": nu_m, "tracing": nu_t, "curve_tracing": [xt, yt], "curve_model": [xm, ym]}))

    # C3b — par paire de billes, corrigé de la distance.
    tb = np.full(mb.shape, np.nan)
    for a in range(len(beads)):
        ra = np.flatnonzero(li == a)
        for b in range(a + 1, len(beads)):
            rb = np.flatnonzero(li == b)
            if len(ra) and len(rb):
                v = t.med[np.ix_(ra, rb)]
                if np.isfinite(v).any():
                    tb[a, b] = np.nanmedian(v)
    k = (bj - bi).astype(np.int64)
    data_oe = oe_vec(tb[bi, bj], k)
    r = spearman(oe_vec(mb[bi, bj], k), data_oe)
    p95, mean = null_p95(lambda v: spearman(v, data_oe),
                         lambda rng: (lambda p: oe_vec(mb[np.ix_(p, p)][bi, bj], k))(rng.permutation(len(beads))))
    res.append(Result("C3b", "Su 2020 chr21 (IMR90)", "ρ corrigé de la distance, par paire de billes",
                      r, "≥ 0,3", verdict(r, 0.3, p95=p95), null_p95=p95,
                      extra={"null_mean": mean, "bead_pairs": int(np.isfinite(data_oe).sum())}))

    # C3c — échelle absolue.
    ratios = {}
    for s0 in (1e6, 5e6):
        w = (sep >= 0.8 * s0) & (sep <= 1.2 * s0)
        mt, mm = np.nanmedian(t.med[iu, ju][w]), np.nanmedian(ml[iu, ju][w])
        ratios[f"{s0 / 1e6:.0f} Mb"] = {"model_um": float(mm), "tracing_um": float(mt), "ratio": float(mm / mt)}
    ok = all(1 / 1.5 <= v["ratio"] <= 1.5 for v in ratios.values())
    worst = max(ratios.values(), key=lambda v: abs(np.log(v["ratio"])))["ratio"]
    res.append(Result("C3c", "Su 2020 chr21 (IMR90)", "rapport des distances médianes modèle / traçage à 1 et 5 Mb",
                      worst, "entre 1/1,5 et 1,5", "accord" if ok else "échec", extra=ratios))
    return res


@dataclass
class Genome:
    chrom: np.ndarray      # (1041,) chromosome de chaque locus
    pos: np.ndarray
    radial: np.ndarray     # médiane sur les cellules de la position radiale normalisée
    lamina: np.ndarray     # médiane de la distance à la lamina (µm), donnée par les auteurs
    med: np.ndarray        # (1041, 1041) µm
    cells: int


def load_genome() -> Genome:
    df = pd.read_csv(locked(manifest(), "su2020.genome").path, sep="\t")
    df.columns = ["z", "x", "y", "coord", "homolog", "cell", "exp", "lamina"]
    chrom, mid = _mid(df["coord"])
    order = sorted(set(zip(chrom, mid)), key=lambda cp: (CHROMS.index(cp[0]) if cp[0] in CHROMS else 99, cp[1]))
    key = {cp: n for n, cp in enumerate(order)}
    li = np.array([key[cp] for cp in zip(chrom, mid)])
    cells, ci = np.unique(df["exp"].astype(str) + ":" + df["cell"].astype(str), return_inverse=True)
    hi = df["homolog"].to_numpy() - 1
    T = np.full((len(cells), len(order), 2, 3), np.nan, dtype=np.float32)
    T[ci, li, hi] = df[["z", "x", "y"]].to_numpy(dtype=np.float32)
    L = np.full((len(cells), len(order), 2), np.nan, dtype=np.float32)
    L[ci, li, hi] = df["lamina"].to_numpy(dtype=np.float32)

    # Position radiale : distance au centroïde des loci détectés de la cellule, sur le 95ᵉ centile.
    flat = T.reshape(len(cells), -1, 3)
    c = np.nanmean(flat, axis=1)
    d = np.linalg.norm(T - c[:, None, None, :], axis=3)
    r95 = np.nanpercentile(d.reshape(len(cells), -1), 95, axis=1)
    rad = np.nanmedian((d / r95[:, None, None]).transpose(1, 0, 2).reshape(len(order), -1), axis=1)
    lam = np.nanmedian(L.transpose(1, 0, 2).reshape(len(order), -1), axis=1) / 1000.0
    med = _block_median(T, np.arange(len(order)), block=3)
    return Genome(np.array([o[0] for o in order]), np.array([o[1] for o in order]), rad, lam, med, len(cells))


def genome_criteria(model: Model, g: Genome) -> tuple[list[Result], dict]:
    res = []
    rb = model.radial()
    idx = np.full((2, len(g.pos)), -1, dtype=np.int64)
    for ch in np.unique(g.chrom):
        sel = np.flatnonzero(g.chrom == ch)
        idx[:, sel] = model.beads_at(ch, g.pos[sel])[:2]
    placed = (idx >= 0).all(axis=0)
    mr = np.where(placed, rb[np.maximum(idx, 0)].mean(axis=0), np.nan)
    groups = [np.flatnonzero(g.chrom == ch) for ch in np.unique(g.chrom)]

    def permute(rng):
        p = np.arange(len(g.pos))
        for s in groups:
            p[s] = rng.permutation(s)
        return p

    r = spearman(mr, g.radial)
    p95, mean = null_p95(lambda v: spearman(v, g.radial), lambda rng: mr[permute(rng)])
    res.append(Result("C4a", "Su 2020, génome (IMR90)", "ρ des positions radiales, par locus", r, "≥ 0,3",
                      verdict(r, 0.3, p95=p95), null_p95=p95,
                      extra={"loci": int(np.isfinite(mr).sum()), "null_mean": mean,
                             "lamina_rho": spearman(mr, -g.lamina)}))

    names = [ch for ch in CHROMS if (g.chrom == ch).any()]
    pt = np.array([np.nanmean(g.radial[g.chrom == ch]) for ch in names])
    pm = np.array([np.nanmean(mr[g.chrom == ch]) for ch in names])
    r = spearman(pm, pt)
    res.append(Result("C4b", "Su 2020, génome (IMR90)", "ρ des positions radiales moyennes, par chromosome",
                      r, "≥ 0,5", "accord" if r >= 0.5 else ("partiel" if r > 0 else "échec"),
                      extra={"chromosomes": names, "tracing": pt.tolist(), "model": pm.tolist()}))

    # C4c — distances entre loci de chromosomes différents.
    mm = np.full(g.med.shape, np.nan, dtype=np.float32)
    rows = np.arange(len(g.pos))
    for b0 in range(0, len(rows), 24):
        r_ = rows[b0:b0 + 24]
        D = pair_distances(model, idx[:, r_], cis_only=False, cols=idx)
        mm[r_] = np.nanmedian(D, axis=0)
    mm[:, ~placed] = np.nan; mm[~placed] = np.nan
    iu, ju = np.triu_indices(len(g.pos), 1)
    inter = g.chrom[iu] != g.chrom[ju]
    iu, ju = iu[inter], ju[inter]
    r = spearman(mm[iu, ju], g.med[iu, ju])
    p95, mean = null_p95(lambda v: spearman(v, g.med[iu, ju]),
                         lambda rng: (lambda p: mm[p][:, p][iu, ju])(permute(rng)))
    res.append(Result("C4c", "Su 2020, génome (IMR90)", "ρ des distances entre loci de chromosomes différents",
                      r, "≥ 0,3", verdict(r, 0.3, p95=p95), null_p95=p95,
                      extra={"pairs": int(len(iu)), "null_mean": mean}))
    figs = {"radial_model": mr.tolist(), "radial_tracing": g.radial.tolist(), "chrom": g.chrom.tolist()}
    return res, figs
