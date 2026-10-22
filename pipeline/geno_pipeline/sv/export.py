"""Ce que la page du variant structural lit : `epha4.variant.json`.

Le locus (gènes, sites CTCF orientés, repères reconstruits et leur source), la famille de
variants, les paramètres d'extrusion et la loi calée — de quoi refaire le 1D en direct dans le
navigateur —, une conformation représentative par allèle simulé en 3D, alignée sur le type
sauvage par les monomères qu'ils partagent, et les couches cliniques du locus (ClinVar P/LP,
GWAS, contrainte gnomAD), par monomère de 2 kb, pour la surimpression.
"""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path

import numpy as np

from ..ensemble.stats import kabsch
from . import locus as L
from .check import PARAMS
from .predict import calibrate
from .run import alleles, load

FORMAT = "geno-variant/1"


def medoid(coords: np.ndarray, sample: int = 120) -> int:
    """La conformation la plus proche de toutes les autres (RMSD après rotation), parmi
    `sample` instantanés pris régulièrement — un tirage typique, pas une moyenne qui n'existe
    pas."""
    pick = np.linspace(0, len(coords) - 1, min(sample, len(coords))).astype(int)
    x = coords[pick] - coords[pick].mean(axis=1, keepdims=True)
    d = np.zeros((len(pick), len(pick)))
    for i in range(len(pick)):
        for j in range(i + 1, len(pick)):
            d[i, j] = d[j, i] = kabsch(x[i], x[j])[1]
    return int(pick[d.sum(axis=1).argmin()])


def _f32(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a, "<f4").tobytes()).decode()


def clinical(loc: L.Locus) -> dict:
    """Par monomère du type sauvage : variants ClinVar pathogènes ou probablement pathogènes,
    associations GWAS ; par gène du locus : pLI et LOEUF de gnomAD."""
    from ..fiche import sources as S

    m = S.manifest()
    n = loc.region.n
    plp = np.zeros(n, np.int32)
    cv, _ = S.clinvar(m, [])
    p = cv.pos[L.CHROM]
    lo, hi = np.searchsorted(p, loc.region.start), np.searchsorted(p, loc.region.end)
    for pos, i in zip(p[lo:hi], cv.idx[L.CHROM][lo:hi]):
        item = cv.items[int(i)]
        sig = item if isinstance(item, str) else item.get("sig")
        if sig == "plp":
            plp[int((pos - loc.region.start) // loc.region.bp_per_bead)] += 1
    gw = S.gwas(m)
    gwas = np.zeros(n, np.int32)
    q = gw.pos.get(L.CHROM, np.empty(0, np.int64))
    lo, hi = np.searchsorted(q, loc.region.start), np.searchsorted(q, loc.region.end)
    for pos in q[lo:hi]:
        gwas[int((pos - loc.region.start) // loc.region.bp_per_bead)] += 1
    genes_iv, _ = S.genes(m)
    gn = S.gnomad(m)
    constraint = []
    for s, e, i in zip(genes_iv.starts[L.CHROM], genes_iv.ends[L.CHROM], genes_iv.idx[L.CHROM]):
        g = genes_iv.items[int(i)]
        if g["type"] != "protein_coding" or not (loc.region.start <= s < loc.region.end):
            continue
        pli, loeuf = gn.get(g["id"], (float("nan"), float("nan")))
        constraint.append({"name": g["name"], "id": g["id"], "start": int(s), "end": int(e), "strand": g["strand"],
                           "pLI": None if np.isnan(pli) else round(pli, 4),
                           "LOEUF": None if np.isnan(loeuf) else round(loeuf, 3)})
    return {"clinvar_plp": plp.tolist(), "gwas": gwas.tolist(), "genes": constraint,
            "sources": {k: m[k].sha256 for k in ("clinvar.vcf", "gwas.catalog", "gnomad.constraint", "gencode.annotation.gtf")}}


def clinical_check(loc: L.Locus, clin: dict, n: int = 200, seed: int = 15) -> dict:
    """Les comptes par monomère, refaits **par un autre chemin** : le balayage ligne à ligne des
    fichiers bruts de la semaine 13 (`fiche.check._scan_points`), sur tous les monomères qui ont
    un compte et `n` autres tirés au hasard. Témoin : les mêmes intervalles décalés d'un monomère."""
    import random

    from ..fiche import sources as S
    from ..fiche.check import _scan_points

    m = S.manifest()
    bp = loc.region.bp_per_bead
    plp, gw = np.asarray(clin["clinvar_plp"]), np.asarray(clin["gwas"])
    ks = sorted(set(np.flatnonzero(plp | gw).tolist()) | set(random.Random(seed).sample(range(loc.region.n), n)))
    out = {"monomers": len(ks)}
    for shift in (0, 1):
        iv = [(L.CHROM, loc.region.start + (k + shift) * bp, loc.region.start + (k + shift + 1) * bp) for k in ks]
        cv = _scan_points(S.locked(m, "clinvar.vcf").path, iv, "clinvar")
        g = _scan_points(S.locked(m, "gwas.catalog").path, iv, "gwas")
        out["shifted" if shift else "equal"] = {
            "clinvar": sum(int(c.get("plp", 0) == plp[k]) for c, k in zip(cv, ks)),
            "gwas": sum(int(x.get("n", 0) == gw[k]) for x, k in zip(g, ks)),
        }
    out["nonzero"] = {"clinvar": int((plp > 0).sum()), "gwas": int((gw > 0).sum())}
    return out


def export(sv_dir: Path, out: Path) -> dict:
    loc = L.build()
    report = sv_dir / "check.json"
    if report.exists():
        law = json.loads(report.read_text())["law"]
    else:
        coords, legs, _ = load(sv_dir / "epha4-wt.npz")
        law = calibrate(coords, legs).__dict__
    al = {}
    ref = None
    ref_mids = None
    for tag, v in alleles(loc).items():
        path = sv_dir / f"epha4-{tag}.npz"
        if not path.exists():
            continue
        coords, _, meta = load(path)
        k = medoid(coords)
        x = coords[k].astype(np.float64)
        mids = np.asarray(meta["mids_hg38"])
        if ref is None:
            ref, ref_mids = x - x.mean(axis=0), mids
            x = ref
        else:
            # Aligné sur le type sauvage par les monomères communs : le fondu de l'un à l'autre
            # ne doit montrer que ce que la délétion change.
            common, ia, ib = np.intersect1d(ref_mids, mids, return_indices=True)
            a, b = ref[ia], x[ib]
            rot, _ = kabsch(a, b)
            x = (x - b.mean(axis=0)) @ rot + a.mean(axis=0)
        al[tag] = {"name": "type sauvage" if v is None else v.name, "kind": None if v is None else v.kind,
                   "variant": None if v is None else {"start": v.start, "end": v.end, "keep": [list(q) for q in v.keep]},
                   "frame": k, "frames": int(len(coords)), "n": int(coords.shape[1]),
                   "sigma_nm": meta["sigma_nm"], "coords": _f32(x), "mids_first": int(mids[0])}
    clin = clinical(loc)
    clin["check"] = clinical_check(loc, clin)
    doc = {
        "format": FORMAT,
        "assembly": "GRCh38", "cell_type": "GM12878",
        "region": {"chrom": loc.region.chrom, "start": loc.region.start, "end": loc.region.end,
                   "bp_per_bead": loc.region.bp_per_bead, "n": loc.region.n},
        "genes": [{"name": g.name, "start": g.start, "end": g.end, "strand": g.strand, "tss": g.tss} for g in loc.genes.values()],
        "ctcf": [{"start": p.start, "end": p.end, "motif": p.motif, "signal": round(p.signal, 1),
                  "occupancy": None if p.motif is None else round(loc.occupancy(p.signal), 4)} for p in loc.peaks],
        "enhancer": list(loc.enhancer), "boundary": list(loc.boundary), "tad": list(loc.tad),
        "promoter_half": L.PROMOTER_HALF,
        "variants": [{"name": v.name, "kind": v.kind, "start": v.start, "end": v.end, "keep": [list(q) for q in v.keep]}
                     for v in loc.variants],
        "params": PARAMS, "law": law,
        "alleles": al,
        "clinical": clin,
        "evidence": {
            "genes": "annotation (GENCODE v50)", "ctcf": "mesuré (ENCODE GM12878, orienté MA0139.1)",
            "tad": "mesuré (Rao 2014, boucle levée hg19 → hg38)",
            "enhancer": "reconstruit depuis le texte de Lupiáñez 2015 (1,66 Mb du promoteur d'EPHA4 ; distance mesurée chez la souris)",
            "variants": "reconstruits depuis le texte (1,75–1,9 Mb ; Table S1 non lue)",
            "structures": "simulées (OpenMM, 8 × 150 conformations par allèle)",
            "prediction": "extrusion 1D + loi calée sur le 3D du type sauvage",
        },
        "source": L.SOURCE,
    }
    raw = json.dumps(doc, ensure_ascii=False, separators=(",", ":")).encode()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(raw)
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(), "alleles": list(al), "clinical_check": clin["check"]}
