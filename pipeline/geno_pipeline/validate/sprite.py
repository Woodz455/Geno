"""C1 et C2 : le modèle contre SPRITE dans GM12878 (Quinodoz et al. 2018), à 1 Mb.

Les cartes SPRITE sont en hg19 : le centre de chaque bin est levé en hg38 par la chaîne Ensembl
verrouillée, et ce qui ne se lève pas sur le même chromosome est compté et écarté. Un bin sans
couverture SPRITE (ligne nulle) est écarté aussi.

Cartes utilisées :
- **intra** : les 23 matrices `human_chrN_1Mb.txt` de GSE114242, telles que publiées (normalisées
  par Hi-Corrector, contacts pondérés 2/n) ;
- **inter** : `human_inter_1Mb_nover2_1000_raw.txt`, pondération 2/n et amas de 1 000 lectures au
  plus. C'est la pondération de l'article ; les comptes bruts suffisent, puisque l'enrichissement
  se rapporte aux marges.
"""

from __future__ import annotations

import gzip
import io
import tarfile
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..ensemble.stats import contact_pairs, slope
from ..fiche.sources import locked, manifest
from ..liftover import Chain
from .common import (CUTOFF, Model, Result, null_p95, pair_distances, spearman, verdict)

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
BIN = 1_000_000
INTER = "human_inter_1Mb_nover2_1000_raw.txt"


def oe_vec(v: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Observé / attendu d'un vecteur de paires : chaque valeur sur la moyenne de son décalage."""
    ok = np.isfinite(v)
    s = np.bincount(k[ok], weights=v[ok], minlength=k.max() + 1)
    c = np.bincount(k[ok], minlength=k.max() + 1)
    mu = np.divide(s, c, out=np.full(len(s), np.nan), where=c > 0)
    return v / mu[k]


@dataclass
class Chrom:
    """Un chromosome : ses bins valides, ses paires, les deux côtés."""

    name: str
    pos38: np.ndarray          # (n,) centres levés en hg38
    i: np.ndarray              # paires (i < j) retenues
    j: np.ndarray
    k: np.ndarray              # décalage en bins hg19
    sprite: np.ndarray         # valeur SPRITE de chaque paire
    med: np.ndarray            # (n, n) distance médiane du modèle
    contact: np.ndarray        # (n, n) probabilité de contact du modèle
    genes: np.ndarray          # (n,) gènes codants à ±500 kb du centre
    matrix: np.ndarray         # (n, n) SPRITE, bins valides
    lost: int = 0
    same_bead: int = 0
    extra: dict = field(default_factory=dict)


def gene_starts() -> dict[str, np.ndarray]:
    """Débuts des gènes codants GENCODE, hg38, par chromosome (pour orienter E1)."""
    path = locked(manifest(), "gencode.annotation.gtf").path
    out: dict[str, list[int]] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 9)
            if p[2] != "gene" or 'gene_type "protein_coding"' not in p[8]:
                continue
            out.setdefault(p[0], []).append(int(p[3]) - 1 if p[6] == "+" else int(p[4]) - 1)
    return {c: np.sort(np.array(v)) for c, v in out.items()}


def load_intra(model: Model, chain: Chain, genes: dict[str, np.ndarray]) -> list[Chrom]:
    m = manifest()
    out = []
    with tarfile.open(locked(m, "sprite.gm12878.intra1mb").path) as tar:
        for chrom in CHROMS:
            raw = tar.extractfile(f"human_{chrom}_1Mb.txt").read().decode()
            M = pd.read_csv(io.StringIO(raw), sep="\t", header=None).to_numpy(dtype=np.float64)
            n = M.shape[0]
            centres19 = (np.arange(n) * BIN + BIN // 2).astype(np.int64)
            cs, ps, _ = chain.points(chrom, centres19)
            lifted = np.array([c == chrom for c in cs]) & (ps >= 0)
            covered = M.sum(axis=1) > 0
            keep = lifted & covered
            idx_all = model.beads_at(chrom, np.where(keep, ps, 0))
            keep &= (idx_all >= 0).all(axis=0)
            sel = np.flatnonzero(keep)
            idx = idx_all[:, sel]
            pos38 = ps[sel]
            D = pair_distances(model, idx, cis_only=False)
            rr = model.radius[np.maximum(idx[0], 0)]
            touch = (D < CUTOFF * (rr[:, None] + rr[None, :])[None]).astype(np.float32)
            touch[np.isnan(D)] = np.nan
            med = np.nanmedian(D, axis=0)
            contact = np.nanmean(touch, axis=0)
            iu, ju = np.triu_indices(len(sel), 1)
            same = idx[0][iu] == idx[0][ju]
            iu, ju = iu[~same], ju[~same]
            Ms = M[np.ix_(sel, sel)]
            g = genes.get(chrom, np.array([], dtype=np.int64))
            gd = np.searchsorted(g, pos38 + BIN // 2) - np.searchsorted(g, pos38 - BIN // 2)
            out.append(Chrom(chrom, pos38, iu, ju, (sel[ju] - sel[iu]).astype(np.int64), Ms[iu, ju],
                             med, contact, gd.astype(float), Ms, lost=int((~lifted).sum()), same_bead=int(same.sum())))
    return out


def e1(matrix: np.ndarray, genes: np.ndarray) -> np.ndarray:
    """Premier vecteur propre de la corrélation de l'observé/attendu, orienté par les gènes."""
    n = matrix.shape[0]
    k = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    v = matrix.astype(float).copy()
    np.fill_diagonal(v, np.nan)
    oe = v / np.array([np.nanmean(v[k == d]) if np.isfinite(v[k == d]).any() else np.nan for d in range(n)])[k]
    np.fill_diagonal(oe, 1.0)
    oe = np.nan_to_num(oe, nan=1.0)
    c = np.corrcoef(oe)
    c = np.nan_to_num(c)
    w, vec = np.linalg.eigh(c)
    e = vec[:, -1]
    if np.corrcoef(e, genes)[0, 1] < 0:
        e = -e
    return e


def criteria(model: Model, chroms: list[Chrom], inter: "Inter") -> list[Result]:
    res: list[Result] = []
    # C1a — brut, contre la distance seule (N1 : la loi moyenne du modèle à chaque décalage).
    prox = np.concatenate([-c.med[c.i, c.j] for c in chroms])
    data = np.concatenate([c.sprite for c in chroms])
    k = np.concatenate([c.k for c in chroms])
    ok = np.isfinite(prox) & np.isfinite(data)
    s = np.bincount(k[ok], weights=prox[ok]); n = np.bincount(k[ok])
    n1 = (s / np.maximum(n, 1))[k]
    r_model, r_n1 = spearman(prox, data), spearman(n1, data)
    res.append(Result("C1a", "SPRITE GM12878, intra, 1 Mb", "ρ brut, proximité du modèle contre SPRITE",
                      r_model, "> ρ(N1)", "accord" if r_model > r_n1 else "échec", reference=r_n1,
                      extra={"pairs": int(ok.sum())}))

    # C1b — corrigé de la distance, contre les permutations à l'intérieur de chaque chromosome.
    data_oe = np.concatenate([oe_vec(c.sprite, c.k) for c in chroms])

    def model_oe(perms=None):
        parts = []
        for t, c in enumerate(chroms):
            med = c.med if perms is None else c.med[np.ix_(perms[t], perms[t])]
            parts.append(oe_vec(1.0 / med[c.i, c.j], c.k))   # proximité = 1 / distance
        return np.concatenate(parts)

    r = spearman(model_oe(), data_oe)
    p95, mean = null_p95(lambda v: spearman(v, data_oe),
                         lambda rng: model_oe([rng.permutation(c.med.shape[0]) for c in chroms]))
    res.append(Result("C1b", "SPRITE GM12878, intra, 1 Mb", "ρ corrigé de la distance", r, "≥ 0,3",
                      verdict(r, 0.3, p95=p95), null_p95=p95, extra={"null_mean": mean}))

    # C1c — compartiments.
    agree = total = 0
    per = {}
    for c in chroms:
        a = e1(c.matrix, c.genes)
        b = e1(1.0 / np.where(c.med > 0, c.med, np.nan), c.genes)
        same = np.sign(a) == np.sign(b)
        agree += int(same.sum()); total += len(a)
        per[c.name] = float(same.mean())
    frac = agree / max(total, 1)
    res.append(Result("C1c", "SPRITE GM12878, intra, 1 Mb", "signe de E1, modèle = SPRITE (part des bins)",
                      frac, "≥ 70 %", "accord" if frac >= 0.7 else ("partiel" if frac > 0.5 else "échec"),
                      reference=0.5, extra={"bins": total, "per_chrom": per}))

    # C1d — pente de P(s), 1–10 Mb.
    sp = np.concatenate([c.sprite for c in chroms]); mc = np.concatenate([c.contact[c.i, c.j] for c in chroms])
    def curve(v):
        okv = np.isfinite(v)
        s = np.bincount(k[okv], weights=v[okv]); n = np.bincount(k[okv])
        return np.arange(len(s)) * BIN, s / np.maximum(n, 1)
    xs, ys = curve(sp); xm, ym = curve(mc)
    a_s, a_m = slope(xs, ys, 1e6, 1e7), slope(xm, ym, 1e6, 1e7)
    res.append(Result("C1d", "SPRITE GM12878, intra, 1 Mb", "pente de P(s), 1–10 Mb : modèle − SPRITE",
                      a_m - a_s, "|écart| ≤ 0,2", "accord" if abs(a_m - a_s) <= 0.2 else "échec",
                      extra={"model": a_m, "sprite": a_s, "curve_sprite": [xs[:40].tolist(), ys[:40].tolist()],
                             "curve_model": [xm[:40].tolist(), ym[:40].tolist()]}))

    # C2 — paires de chromosomes.
    e_model, e_data = inter.model, inter.data
    iu = np.triu_indices(len(CHROMS), 1)
    r = spearman(e_model[iu], e_data[iu])
    # N2 ne s'applique pas ici : permuter les billes à l'intérieur d'un chromosome ne change aucune
    # somme par paire de chromosomes. La référence permute les chromosomes entre eux.
    p95, mean = null_p95(lambda v: spearman(v, e_data[iu]),
                         lambda rng: (lambda p: e_model[np.ix_(p, p)][iu])(rng.permutation(len(CHROMS))))
    res.append(Result("C2", "SPRITE GM12878, inter, 1 Mb", "ρ des enrichissements par paire de chromosomes (253)",
                      r, "≥ 0,3", verdict(r, 0.3, p95=p95), null_p95=p95,
                      extra={"null_mean": mean, "null": "chromosomes permutés entre eux (voir le module)"}))
    return res


@dataclass
class Inter:
    data: np.ndarray    # (23, 23) enrichissement SPRITE
    model: np.ndarray   # (23, 23) enrichissement du modèle


def enrichment(M: np.ndarray) -> np.ndarray:
    M = M.astype(float).copy()
    np.fill_diagonal(M, 0.0)
    R = M.sum(axis=1)
    T = R.sum()
    E = np.divide(M * T, np.outer(R, R), out=np.zeros_like(M), where=np.outer(R, R) > 0)
    np.fill_diagonal(E, np.nan)
    return E


def load_inter(model: Model, sizes: list[int]) -> Inter:
    with tarfile.open(locked(manifest(), "sprite.gm12878.inter1mb").path) as tar:
        M = pd.read_csv(tar.extractfile(INTER), sep="\t", header=None, dtype=np.float64).to_numpy()
    edges = np.concatenate([[0], np.cumsum(sizes)])
    if M.shape[0] != edges[-1]:
        raise ValueError(f"carte inter : {M.shape[0]} bins, attendu {edges[-1]}")
    blocks = np.add.reduceat(np.add.reduceat(M, edges[:-1], axis=0), edges[:-1], axis=1)

    chrom_of_bead = np.array([CHROMS.index(model.chrom_of_copy(c)) for c in model.copy_id])
    C = np.zeros((len(CHROMS), len(CHROMS)))
    for x in model.coords:
        p = contact_pairs(x.astype(np.float64), model.radius, CUTOFF)
        a, b = chrom_of_bead[p[:, 0]], chrom_of_bead[p[:, 1]]
        inter = a != b
        np.add.at(C, (a[inter], b[inter]), 1)
        np.add.at(C, (b[inter], a[inter]), 1)
    return Inter(enrichment(blocks), enrichment(C))
