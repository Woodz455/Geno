"""Compose les fiches de toutes les billes d'un `.g3d`, et les écrit dans un fichier à part.

**Ce qu'est une fiche.** Pour une bille — ou un monomère, ou une base — tout ce que l'on
sait de l'intervalle génomique qu'elle porte, rangé dans l'ordre de la feuille de route :
coordonnées, structure, gènes, régulation, clinique, séquence. Chaque champ dit **d'où il
vient** (`s`, un identifiant du manifeste) et **ce qu'il vaut** (`e`) :

    mesure      mesuré dans GM12878
    autre       mesuré, mais dans une autre lignée ou un autre tissu (`cell` dit lequel)
    simule      produit par le modèle — l'ensemble des semaines 7–8
    reference   séquence ou assemblage de référence
    annotation  annotation de référence (GENCODE, ClinVar, GWAS, gnomAD…)

Un champ sans valeur porte `absent`, une phrase qui dit pourquoi. Aucun champ n'est muet :
c'est ce que mesure `check`, et ce que veut dire « fiche complète ».

**Deux parties.** Ce qui ne dépend que de l'intervalle (les annotations) est partagé par
les homologues : un enregistrement par intervalle distinct, rangés par chromosome et par
blocs de 64. Ce qui dépend de la copie (E1 et position radiale du modèle) tient dans une
colonne par bille. Le navigateur assemble les deux.

**Le fichier.** Même principe que le `.g3d` (`export/g3d.py`) : préambule, en-tête JSON
compressé, puis des blocs deflate qu'on lit par requête Range, chacun avec son empreinte.
Il est à part pour que le premier rendu de la semaine 9 n'en paie rien.
"""

from __future__ import annotations

import hashlib
import json
import math
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..export.g3d import _deflate, _inflate, digest
from ..export.g3d import read as read_g3d
from ..sequence import TwoBit
from . import sources as S

MAGIC = b"\x89FCH\r\n\x1a\n"
VERSION = 1
PREAMBLE = 64
BLOCK = 64
TOP = 5             # au-delà, un compte et « +k » (les niveaux fins, eux, listent tout)

EVIDENCE = {
    "mesure": "mesuré dans GM12878",
    "autre": "mesuré dans une autre lignée ou un autre tissu",
    "simule": "simulé — produit par le modèle",
    "reference": "séquence ou assemblage de référence",
    "annotation": "annotation de référence",
}

SECTIONS = {
    "coordonnees": ["hg38", "t2t", "cytobande"],
    "structure": ["compartiment", "sous_compartiment", "tad", "boucles", "radial",
                  "lad_modele", "lad", "speckles", "nad", "replication"],
    "genes": ["genes", "expression", "contrainte"],
    "regulation": ["ccre", "chromhmm", "ctcf"],
    "clinique": ["clinvar", "gwas"],
    "sequence": ["composition", "bases"],
}
FINE_ONLY = {"structure": ["ctcf_plante", "boucles_modele"]}

LABELS = {
    "hg38": "Position hg38", "t2t": "Position T2T-CHM13", "cytobande": "Cytobande",
    "compartiment": "Compartiment (modèle)", "sous_compartiment": "Sous-compartiment",
    "tad": "TAD", "boucles": "Ancres de boucles", "radial": "Position radiale",
    "lad_modele": "LAD (modèle)", "lad": "Lamina (LMNB1)", "speckles": "Speckles (SON)",
    "nad": "Nucléole (NAD)", "replication": "Timing de réplication",
    "ctcf_plante": "CTCF du modèle", "boucles_modele": "Boucles du modèle",
    "genes": "Gènes", "expression": "Expression", "contrainte": "Contrainte LoF",
    "ccre": "cCRE", "chromhmm": "ChromHMM", "ctcf": "CTCF (ChIP-seq)",
    "clinvar": "ClinVar", "gwas": "GWAS", "composition": "Composition", "bases": "Bases",
}


def F(e: str, s: str, v=None, *, absent: str | None = None, **kw) -> dict:
    d = {"e": e, "s": s, **{k: x for k, x in kw.items() if x is not None}}
    if absent is not None:
        d["absent"] = absent
    else:
        d["v"] = v
    return d


def _r(x: float, k: int = 3):
    return None if x is None or not math.isfinite(x) else round(float(x), k)


# ---------------------------------------------------------------------------
# Le contexte : toutes les sources, lues une fois
# ---------------------------------------------------------------------------

@dataclass
class Context:
    m: dict
    genome: TwoBit
    hs1: object
    cyto: S.Intervals
    genes: S.Intervals
    gencode: dict
    pc_starts: dict
    gtex: dict
    tissues: list
    gnomad: dict
    ccre: S.Intervals
    chromhmm: S.Intervals
    ctcf: S.Intervals
    subcomp: S.Intervals
    domains: S.Intervals
    boundaries: dict
    loops: S.Intervals
    repli: dict
    repli_terciles: tuple[float, float]
    lad: S.BigWig
    son: S.BigWig
    nad_nuc: dict
    nad_ctl: dict
    clinvar: S.Points
    clinvar_release: str
    gwas: S.Points
    lift_report: dict


def context(fine_regions: list[tuple[str, int, int]], log=print) -> Context:
    m = S.manifest()
    log("  sources : coordonnées et chaînes")
    genome = TwoBit(S.locked(m, "hg38.2bit").path)
    hs1 = S.chain_hs1(m)
    c19 = S.chain_hg19(m)
    log("  sources : gènes (GENCODE, GTEx, gnomAD)")
    genes, gencode = S.genes(m)
    pc = {}
    for chrom in genes.starts:
        items = [genes.items[i] for i in genes.idx[chrom]]
        rows = sorted(((g["start"], g["end"], g) for g in items if g["type"] == "protein_coding"),
                      key=lambda r: (r[0], r[1]))
        pc[chrom] = (np.array([r[0] for r in rows]), np.array([r[1] for r in rows]), [r[2] for r in rows])
    gt, tissues = S.gtex(m)
    log("  sources : régulation (cCRE, ChromHMM, CTCF + motif)")
    ctcf = S.ctcf(m, genome)
    log("  sources : Rao 2014 (levé depuis hg19)")
    sub, sub_i = S.subcompartments(m, c19)
    dom, dom_i = S.domains(m, c19)
    bounds = {}
    for chrom in dom.starts:
        bounds[chrom] = np.unique(np.concatenate([dom.starts[chrom], dom.ends[chrom]]))
    lp, lp_i = S.loops(m, c19)
    log("  sources : repères nucléaires")
    repli, rep_i = S.replication(m, c19)
    allv = np.concatenate([v for _, v in repli.values()])
    t1, t2 = np.quantile(allv, [1 / 3, 2 / 3])
    nad_n = S.bedgraph(S.locked(m, "nads.hffc6.nucleolus").path)
    nad_c = S.bedgraph(S.locked(m, "nads.hffc6.control").path)
    log("  sources : ClinVar, GWAS")
    cv, cv_i = S.clinvar(m, fine_regions)
    gw = S.gwas(m)
    return Context(
        m=m, genome=genome, hs1=hs1, cyto=S.cytobands(m), genes=genes, gencode=gencode,
        pc_starts=pc, gtex=gt, tissues=tissues, gnomad=S.gnomad(m), ccre=S.ccres(m),
        chromhmm=S.chromhmm(m), ctcf=ctcf, subcomp=sub, domains=dom, boundaries=bounds,
        loops=lp, repli=repli, repli_terciles=(float(t1), float(t2)),
        lad=S.BigWig(S.locked(m, "lads.damid.lmnb1").path), son=S.BigWig(S.locked(m, "tsaseq.son").path),
        nad_nuc=nad_n, nad_ctl=nad_c, clinvar=cv, clinvar_release=cv_i["release"], gwas=gw,
        lift_report={"rao2014.gm12878.subcompartments": sub_i, "rao2014.gm12878.domains": dom_i,
                     "rao2014.gm12878.loops": lp_i, "repliseq.gm12878": rep_i},
    )


# ---------------------------------------------------------------------------
# Un intervalle → son enregistrement
# ---------------------------------------------------------------------------

class Tables:
    """Ce que plusieurs fiches partagent, écrit une fois : gènes (avec expression et
    contrainte), maladies, traits GWAS. Une fiche n'en garde que l'indice. Mesuré : les
    gènes recopiés dans chaque fiche faisaient à eux seuls 3,2 Mo sur 12 (aperçu + noyau)."""

    def __init__(self, ctx: "Context"):
        self.ctx = ctx
        self.gene_ix: dict[str, int] = {}
        self.genes: list[dict] = []
        # Indices fixés d'avance, par chromosome puis par début : la table s'écrit en un bloc
        # par chromosome, et une fiche n'a besoin que de celui de son chromosome.
        for chrom in sorted(ctx.genes.starts, key=_chrom_key):
            for i in ctx.genes.idx[chrom]:
                g = ctx.genes.items[i]
                if g["id"] not in self.gene_ix:
                    self._add({**g, "chrom": chrom})
        self.dis_ix: dict[str, int] = {}
        self.diseases: list[list] = []
        self.trait_ix: dict[str, int] = {}
        self.traits: list[list] = []

    def _add(self, g: dict) -> int:
        k = self.gene_ix[g["id"]] = len(self.genes)
        x = self.ctx.gtex.get(g["id"])
        c = self.ctx.gnomad.get(g["id"])
        self.genes.append({**{f: g[f] for f in ("id", "name", "type", "strand", "start", "end", "mane", "chrom")},
                           "gtex": [_r(x[0], 2), x[1], _r(x[2], 2)] if x else None,
                           "gnomad": [_r(c[0]), _r(c[1])] if c else None})
        return k

    def gene(self, g: dict) -> int:
        k = self.gene_ix.get(g["id"])
        return self._add(g) if k is None else k

    def blobs(self) -> dict[str, list[tuple[dict, object]]]:
        """Les tables à écrire : gènes en colonnes par chromosome, maladies et traits par
        tranches de 4 096. Rend, par table, des couples (méta, contenu JSON)."""
        out: dict[str, list] = {"genes": [], "diseases": [], "traits": []}
        k = 0
        while k < len(self.genes):
            chrom = self.genes[k]["chrom"]
            j = k
            while j < len(self.genes) and self.genes[j]["chrom"] == chrom:
                j += 1
            part = self.genes[k:j]
            cols = {f: [g[f] for g in part] for f in ("id", "name", "type", "start", "end", "gtex", "gnomad")}
            cols["strand"] = "".join(g["strand"] for g in part)
            cols["mane"] = "".join("1" if g["mane"] else "0" for g in part)
            out["genes"].append(({"chrom": chrom, "first": k, "n": j - k}, cols))
            k = j
        for name, rows in (("diseases", self.diseases), ("traits", self.traits)):
            for k in range(0, len(rows), 4096):
                out[name].append(({"first": k, "n": len(rows[k:k + 4096])}, rows[k:k + 4096]))
        return out

    def disease(self, d: list) -> int:
        ids, name = d[0], d[-1]
        key = "|".join(ids) or name
        k = self.dis_ix.get(key)
        if k is None:
            k = self.dis_ix[key] = len(self.diseases)
            self.diseases.append([ids, name])
        return k

    def trait(self, h: dict) -> int:
        key = h["trait"] + "\t" + h["uri"]
        k = self.trait_ix.get(key)
        if k is None:
            k = self.trait_ix[key] = len(self.traits)
            self.traits.append([h["trait"], h["mapped"], h["uri"]])
        return k

    def variant(self, it: dict, pos: int, fine: bool) -> dict:
        # Hors des régions fines, la fiche ne descend pas à la base : allèles et type de
        # variant y sont déjà dans le HGVS, on ne les répète pas.
        drop = {"diseases"} if fine else {"diseases", "ref", "alt", "kind"}
        return {**{k: v for k, v in it.items() if k not in drop}, "pos": pos,
                "diseases": [self.disease(d) for d in it.get("diseases", [])]}

    def hit(self, h: dict, pos: int) -> dict:
        return {"rs": h["rs"], "pos": pos, "t": self.trait(h), "p": h["p"], "study": h["study"]}


def _nearest(ctx: Context, tb: Tables, chrom: str, s: int, e: int) -> list[dict]:
    if chrom not in ctx.pc_starts:
        return []
    st, en, gs = ctx.pc_starts[chrom]
    out = []
    left = np.nonzero(en <= s)[0]
    if left.size:
        k = left[np.argmax(en[left])]
        out.append({"g": tb.gene({**gs[k], "chrom": chrom}), "distance": int(s - en[k]), "side": "amont"})
    right = np.nonzero(st >= e)[0]
    if right.size:
        k = right[np.argmin(st[right])]
        out.append({"g": tb.gene({**gs[k], "chrom": chrom}), "distance": int(st[k] - e), "side": "aval"})
    return out


def record(ctx: Context, tb: Tables, chrom: str, s: int, e: int, fine: bool) -> dict:
    """L'enregistrement d'un intervalle `[s, e)` de hg38 (0-based, demi-ouvert). Gènes,
    maladies et traits y sont des indices dans `tb` ; `expand` les remet en place."""
    rec: dict[str, dict] = {}
    c = (s + e) // 2

    # --- coordonnées ----------------------------------------------------------
    rec["hg38"] = F("reference", "hg38.chrom.sizes", {"chrom": chrom, "start": s, "end": e})
    iv = ctx.hs1.interval(chrom, s, e)
    if iv is None:
        rec["t2t"] = F("reference", "chain.hg38.to.hs1",
                       absent="aucun bloc de la chaîne hg38 → T2T dans cet intervalle")
    else:
        rec["t2t"] = F("reference", "chain.hg38.to.hs1",
                       {"chrom": iv.chrom, "start": iv.start, "end": iv.end, "strand": iv.strand,
                        "mapped": _r(iv.mapped), "split": iv.split})
    bands = ctx.cyto.get(chrom, s, e)
    rec["cytobande"] = (F("reference", "hg38.cytoband", bands) if bands else
                        F("reference", "hg38.cytoband", absent="hors des bandes publiées par UCSC"))

    # --- structure mesurée ----------------------------------------------------
    sc = {}
    for name, bp in ctx.subcomp.overlap_bp(chrom, s, e):
        sc[name] = sc.get(name, 0) + bp
    if sc:
        tot = e - s
        rec["sous_compartiment"] = F("mesure", "rao2014.gm12878.subcompartments",
                                     [{"name": k, "frac": _r(v / tot)} for k, v in
                                      sorted(sc.items(), key=lambda kv: -kv[1])], lifted="hg19 → hg38")
    else:
        rec["sous_compartiment"] = F("mesure", "rao2014.gm12878.subcompartments", lifted="hg19 → hg38",
                                     absent="non appelé par Rao 2014 ici (NA), ou non levé depuis hg19")

    if chrom in ctx.boundaries:
        b = ctx.boundaries[chrom]
        k = int(np.searchsorted(b, c))
        cand = [b[j] for j in (k - 1, k) if 0 <= j < len(b)]
        nb = min(cand, key=lambda x: abs(int(x) - c))
        doms = ctx.domains.get(chrom, s, e)
        inside = [d for d in ctx.domains.get(chrom, c, c + 1)]
        smallest = min(inside, key=lambda d: d[1] - d[0]) if inside else None
        rec["tad"] = F("mesure", "rao2014.gm12878.domains",
                       {"n": len(doms), "containing": list(smallest) if smallest else None,
                        "nearest_boundary": int(nb), "boundary_distance": int(abs(int(nb) - c))},
                       lifted="hg19 → hg38")
    else:
        rec["tad"] = F("mesure", "rao2014.gm12878.domains", absent="aucun domaine Arrowhead sur ce chromosome")

    lps = sorted(ctx.loops.get(chrom, s, e), key=lambda x: -x["observed"])
    rec["boucles"] = F("mesure", "rao2014.gm12878.loops", {"n": len(lps), "anchors": lps[:TOP]},
                       lifted="hg19 → hg38")

    lad_mean, lad_frac = ctx.lad.summary(chrom, s, e, above=0.0)
    rec["lad"] = (F("autre", "lads.damid.lmnb1", {"log2": _r(lad_mean), "frac": _r(lad_frac)}, cell="K562")
                  if math.isfinite(lad_mean) else
                  F("autre", "lads.damid.lmnb1", cell="K562", absent="pas de signal pA-DamID ici "
                    "(trou d'assemblage, ou chrY absent de K562)"))
    son_mean, _ = ctx.son.summary(chrom, s, e)
    rec["speckles"] = (F("autre", "tsaseq.son", {"score": _r(son_mean)}, cell="K562")
                       if math.isfinite(son_mean) else
                       F("autre", "tsaseq.son", cell="K562", absent="pas de signal TSA-seq ici"))
    nn, nc = S.bedgraph_sum(ctx.nad_nuc, chrom, s, e), S.bedgraph_sum(ctx.nad_ctl, chrom, s, e)
    if math.isfinite(nn) and math.isfinite(nc) and nn > 0 and nc > 0:
        rec["nad"] = F("autre", "nads.hffc6.nucleolus", {"log2": _r(math.log2(nn / nc))}, cell="HFFc6",
                       control="nads.hffc6.control")
    else:
        rec["nad"] = F("autre", "nads.hffc6.nucleolus", cell="HFFc6",
                       absent="pas assez de lectures NAD-seq ici pour un rapport")
    if chrom in ctx.repli:
        pos, val = ctx.repli[chrom]
        half = max((e - s) // 2, S.REPLI_BIN)
        lo, hi = np.searchsorted(pos, c - half), np.searchsorted(pos, c + half)
        if hi > lo:
            v = float(val[lo:hi].mean())
            t1, t2 = ctx.repli_terciles
            rec["replication"] = F("mesure", "repliseq.gm12878", {
                "score": _r(v, 1), "class": "précoce" if v >= t2 else "moyen" if v >= t1 else "tardif"},
                lifted="hg19 → hg38", phases=[k for k, _ in S.REPLI_PHASES])
        else:
            rec["replication"] = F("mesure", "repliseq.gm12878", absent="aucun bin Repli-seq levé ici")
    else:
        rec["replication"] = F("mesure", "repliseq.gm12878",
                               absent="chromosome absent du Repli-seq UW (chrY : GM12878 est femelle)")

    # --- gènes ----------------------------------------------------------------
    gs = sorted(ctx.genes.get(chrom, s, e), key=lambda g: (g["start"], g["name"]))
    src = ctx.gencode["release"]
    if gs:
        rec["genes"] = F("annotation", "gencode.annotation.gtf", {"genes": [tb.gene({**g, "chrom": chrom}) for g in gs]},
                         release=src)
    else:
        rec["genes"] = F("annotation", "gencode.annotation.gtf",
                         {"genes": [], "nearest": _nearest(ctx, tb, chrom, s, e)}, release=src)
    shown = rec["genes"]["v"]["genes"] or [n["g"] for n in rec["genes"]["v"]["nearest"]]
    expr = [k for k in shown if tb.genes[k]["gtex"]]
    rec["expression"] = (F("autre", "gtex.median.tpm", expr, cell="lymphocytes EBV (GTEx v10)") if expr else
                         F("autre", "gtex.median.tpm", cell="lymphocytes EBV (GTEx v10)",
                           absent="aucun de ces gènes n'est dans GTEx (petits ARN, pseudogènes…)"
                           if shown else "aucun gène dans ou près de l'intervalle"))
    cons = [k for k in shown if tb.genes[k]["gnomad"]]
    rec["contrainte"] = (F("annotation", "gnomad.constraint", cons) if cons else
                         F("annotation", "gnomad.constraint",
                           absent="gnomAD ne mesure la contrainte que des gènes codants"))

    # --- régulation -----------------------------------------------------------
    cc = ctx.ccre.get(chrom, s, e)
    counts: dict[str, int] = {}
    low = 0
    for _, cls in cc:
        if cls == "Low-DNase":
            low += 1
        else:
            counts[cls] = counts.get(cls, 0) + 1
    v = {"n": sum(counts.values()), "classes": counts, "inactive": low}
    if fine:
        v["list"] = []
        k = ctx.ccre.query(chrom, s, e)
        if k.size:
            st, en = ctx.ccre.starts[chrom][k], ctx.ccre.ends[chrom][k]
            ids = ctx.ccre.idx[chrom][k]
            v["list"] = [{"id": ctx.ccre.items[i][0], "class": ctx.ccre.items[i][1], "start": int(a), "end": int(b)}
                         for i, a, b in zip(ids, st, en) if ctx.ccre.items[i][1] != "Low-DNase"]
    rec["ccre"] = F("mesure", "encode.ccres", v)

    hm: dict[str, int] = {}
    for st_, bp in ctx.chromhmm.overlap_bp(chrom, s, e):
        hm[st_] = hm.get(st_, 0) + bp
    if hm:
        tot = sum(hm.values())
        v = {"states": [{"state": k, "frac": _r(b / tot)} for k, b in sorted(hm.items(), key=lambda kv: -kv[1])]}
        if fine:
            k = ctx.chromhmm.query(chrom, s, e)
            v["segments"] = [[int(a), int(b), ctx.chromhmm.items[i]] for i, a, b in
                             zip(ctx.chromhmm.idx[chrom][k], ctx.chromhmm.starts[chrom][k], ctx.chromhmm.ends[chrom][k])]
        rec["chromhmm"] = F("mesure", "chromhmm.gm12878", v, model="18 états")
    else:
        rec["chromhmm"] = F("mesure", "chromhmm.gm12878", absent="non segmenté ici (trou d'assemblage)")

    k = ctx.ctcf.query(chrom, s, e)
    peaks = [{**ctx.ctcf.items[i], "start": int(a), "end": int(b)} for i, a, b in
             zip(ctx.ctcf.idx[chrom][k], ctx.ctcf.starts[chrom][k], ctx.ctcf.ends[chrom][k])] if k.size else []
    v = {"n": len(peaks), "plus": sum(p["motif"] == "+" for p in peaks),
         "minus": sum(p["motif"] == "-" for p in peaks), "none": sum(p["motif"] is None for p in peaks)}
    if fine:
        v["list"] = peaks
    rec["ctcf"] = F("mesure", "ctcf.peaks.gm12878", v, motif="ctcf.motif.orientation")

    # --- clinique -------------------------------------------------------------
    ix = ctx.clinvar.range(chrom, s, e)
    pos = ctx.clinvar.pos.get(chrom)
    cnt = {k: 0 for k in S.SIG}
    plp = []
    allv = []
    if ix.size:
        p0 = int(np.searchsorted(pos, s))
        for off, i in enumerate(ix):
            it = ctx.clinvar.items[i]
            sig = it if isinstance(it, str) else it["sig"]
            cnt[sig] += 1
            if not isinstance(it, str):
                d = tb.variant(it, int(pos[p0 + off]), fine)
                if sig == "plp":
                    plp.append(d)
                if fine:
                    allv.append(d)
    plp.sort(key=lambda d: (-d["stars"], d["pos"]))
    # Aux niveaux fins, la liste complète suffit : `expand` en tire les P/LP.
    v = {"counts": cnt, "plp_total": len(plp)}
    if fine:
        v["variants"] = allv
    else:
        v["plp"] = plp[:TOP]
    rec["clinvar"] = F("annotation", "clinvar.vcf", v, release=ctx.clinvar_release,
                       diseases="mondo.obo + hpo.obo ; OMIM : identifiants seulement")

    ix = ctx.gwas.range(chrom, s, e)
    gpos = ctx.gwas.pos.get(chrom)
    hits = []
    if ix.size:
        p0 = int(np.searchsorted(gpos, s))
        hits = sorted(((ctx.gwas.items[i], int(gpos[p0 + off])) for off, i in enumerate(ix)), key=lambda h: h[0]["p"])
    top = hits if fine else hits[:TOP]
    rec["gwas"] = F("annotation", "gwas.catalog", {"n": len(hits), "top": [tb.hit(h, p) for h, p in top]})

    # --- séquence -------------------------------------------------------------
    seq = ctx.genome.fetch(chrom, s, e)
    n_n = seq.count("N")
    acgt = len(seq) - n_n
    gc = (seq.count("G") + seq.count("C")) / acgt if acgt else math.nan
    rec["composition"] = F("reference", "hg38.2bit", {"gc": _r(gc), "n": _r(n_n / max(1, len(seq)))})
    rec["bases"] = (F("reference", "hg38.2bit", {"zoom": True}) if fine else
                    F("reference", "hg38.2bit", absent="les bases se lisent au zoom maximal, dans les "
                      "régions fines (chr7:4–8 Mb, chr17:5,67–9,67 Mb) — pas à cette échelle"))
    return rec


# ---------------------------------------------------------------------------
# Ce qui dépend de la copie : le modèle
# ---------------------------------------------------------------------------

def model_fields(lv: dict, ens: dict | None, parent: dict | None, radius_nm: float) -> list[dict]:
    """Par bille stockée : E1, position radiale (structure affichée, et ensemble si on l'a),
    fraction LAD du modèle — et, aux niveaux fins, les CTCF plantés et les boucles simulées."""
    meta = lv["meta"]
    n = len(lv["start"])
    T = np.asarray(meta["frame"]["transform"], dtype=np.float64).reshape(4, 4)
    p = lv["positions_nm"]
    world = p @ T[:3, :3].T + T[:3, 3]
    radial = np.linalg.norm(world, axis=1) / radius_nm
    copies = meta["copies"]
    e1 = lv["attributes"].get("compartment")
    fine = meta["bp_per_bead"] < 100_000
    out = []
    links = lv.get("links")
    loop_of: dict[int, list] = {}
    if links is not None:
        for i, j, k, w in zip(links["i"], links["j"], links["kind"], links["weight"]):
            for a, b in ((int(i), int(j)), (int(j), int(i))):
                loop_of.setdefault(a, []).append({"partner": b + 1, "kind": links["kinds"][int(k)],
                                                  "weight": _r(w)})
    for b in range(n):
        label = copies[int(lv["copy"][b])]
        f: dict[str, dict] = {}
        src_e1 = "g3d:compartment"
        val = e1[b] if e1 is not None else math.nan
        if fine and parent is not None:
            val = parent["e1"](label, int(lv["start"][b]))
            src_e1 = "g3d:compartment (bille de 750 kb parente)"
        if math.isfinite(val):
            f["compartiment"] = F("simule", src_e1, {"e1": _r(val), "ab": "A" if val > 0 else "B"})
        else:
            f["compartiment"] = F("simule", src_e1, absent="X inactif (pas de compartiments A/B) "
                                  "ou bin non appelé" if label.endswith(":Xi") else "bin non appelé par le modèle")
        rv = {"shown": _r(radial[b])}
        if ens is not None and b in ens["radial"]:
            mu, sd = ens["radial"][b]
            rv.update({"mean": _r(mu), "sd": _r(sd), "over": ens["over"]})
        f["radial"] = F("simule", "g3d:positions", rv,
                        note=None if ens is not None else
                        "conformations fines non posées une à une dans le noyau : seule la structure affichée l'est")
        lad = ens["lad"].get(b) if ens is not None else (parent["lad"](label, int(lv["start"][b])) if parent else None)
        f["lad_modele"] = (F("simule", "ensemble:beads.lad", {"frac": _r(lad)}, track="synthétique")
                           if lad is not None and math.isfinite(lad) else
                           F("simule", "ensemble:beads.lad", absent="pas de piste LAD pour cette bille"))
        if fine:
            fw = lv["attributes"].get("ctcf_forward")
            rv_ = lv["attributes"].get("ctcf_reverse")
            sites = {}
            if fw is not None and math.isfinite(fw[b]):
                sites["+"] = _r(fw[b])
            if rv_ is not None and math.isfinite(rv_[b]):
                sites["-"] = _r(rv_[b])
            f["ctcf_plante"] = F("simule", "g3d:ctcf_forward/ctcf_reverse", {"sites": sites},
                                 note="sites plantés (semaine 8), pas mesurés")
            f["boucles_modele"] = F("simule", "g3d:links", {"loops": loop_of.get(b, [])})
        out.append(f)
    return out


def ensemble_model(path: Path, noyau: dict, apercu: dict | None):
    """Rayon nucléaire, et pour les billes du noyau et de l'aperçu : position radiale
    moyenne et écart-type sur l'ensemble, fraction LAD synthétique."""
    from ..ensemble import read as read_ensemble

    e = read_ensemble(path)
    b = e.beads
    R = float(b.nuclear_radius) * 1000.0
    radial = e.radial()
    key = {(b.labels[int(c)], int(s)): i for i, (c, s) in enumerate(zip(b.copy_id, b.start))}
    nmeta = noyau["meta"]
    out_n = {"radial": {}, "lad": {}, "over": f"{e.n_structures} structures"}
    for k in range(len(noyau["start"])):
        i = key.get((nmeta["copies"][int(noyau["copy"][k])], int(noyau["start"][k])))
        if i is None:
            continue
        out_n["radial"][k] = (float(radial[:, i].mean()), float(radial[:, i].std()))
        out_n["lad"][k] = float(b.lad[i])
    out_a = None
    if apercu is not None:
        out_a = {"radial": {}, "lad": {}, "over": f"{e.n_structures} structures"}
        amet = apercu["meta"]
        by_copy: dict[str, list[int]] = {}
        for i, (c, s) in enumerate(zip(b.copy_id, b.start)):
            by_copy.setdefault(b.labels[int(c)], []).append(i)
        for k in range(len(apercu["start"])):
            label = amet["copies"][int(apercu["copy"][k])]
            s, t = int(apercu["start"][k]), int(apercu["end"][k])
            mem = [i for i in by_copy.get(label, []) if b.start[i] < t and b.end[i] > s]
            if not mem:
                continue
            cen = e.coords[:, mem, :].mean(axis=1)
            r = np.linalg.norm(cen, axis=1) / b.nuclear_radius
            out_a["radial"][k] = (float(r.mean()), float(r.std()))
            out_a["lad"][k] = float(np.mean(b.lad[mem]))
    return R, out_n, out_a


def parent_lookup(noyau: dict, ens_n: dict):
    """Pour un monomère fin : l'E1 et la fraction LAD de la bille de 750 kb qui le contient."""
    meta = noyau["meta"]
    e1 = noyau["attributes"].get("compartment")
    by: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for c, label in enumerate(meta["copies"]):
        k = np.nonzero(noyau["copy"] == c)[0]
        o = np.argsort(noyau["start"][k])
        by[label] = (noyau["start"][k][o], noyau["end"][k][o], k[o])

    def find(label, pos):
        st, en, k = by[label]
        j = int(np.searchsorted(st, pos, side="right")) - 1
        return int(k[j]) if j >= 0 and en[j] > pos else None

    return {
        "e1": lambda label, pos: (float(e1[find(label, pos)]) if e1 is not None and find(label, pos) is not None
                                  else math.nan),
        "lad": lambda label, pos: ens_n["lad"].get(find(label, pos)) if find(label, pos) is not None else None,
    }


# ---------------------------------------------------------------------------
# Le fichier
# ---------------------------------------------------------------------------

def _blob(obj) -> bytes:
    return _deflate(json.dumps(obj, ensure_ascii=False, separators=(",", ":")).encode())


def _chrom_key(chrom: str):
    c = chrom.removeprefix("chr")
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def build(g3d_path: Path, out: Path, ensemble_path: Path | None, log=print) -> dict:
    """Toutes les fiches de toutes les billes du `.g3d`. Rend le bilan (tailles, temps)."""
    import time

    t0 = time.time()
    g = read_g3d(g3d_path)
    lvs = g.levels
    fine_names = [k for k, lv in lvs.items() if lv["meta"]["bp_per_bead"] < 100_000]
    fine_regions = []
    for k in fine_names:
        lv = lvs[k]
        chrom = lv["meta"]["copies"][int(lv["copy"][0])].split(":")[0]
        fine_regions.append((chrom, int(lv["start"].min()), int(lv["end"].max())))
    ctx = context(fine_regions, log=log)
    log(f"  sources lues en {time.time() - t0:.0f} s")

    R, ens_n, ens_a = (5000.0, None, None)
    if ensemble_path and Path(ensemble_path).exists() and "noyau" in lvs:
        R, ens_n, ens_a = ensemble_model(ensemble_path, lvs["noyau"], lvs.get("apercu"))
    parent = parent_lookup(lvs["noyau"], ens_n) if ens_n is not None else None

    blobs: list[bytes] = []
    offset = 0

    def put(b: bytes) -> dict:
        nonlocal offset
        spec = {"offset": offset, "length": len(b), "sha256": digest(b)}
        blobs.append(b)
        offset += len(b)
        return spec

    tb = Tables(ctx)
    levels_meta = []
    report = {"levels": {}}
    for name, lv in lvs.items():
        t1 = time.time()
        fine = name in fine_names
        chroms = [lv["meta"]["copies"][int(c)].split(":")[0] for c in lv["copy"]]
        keys = sorted({(ch, int(s), int(e)) for ch, s, e in zip(chroms, lv["start"], lv["end"])},
                      key=lambda k: (_chrom_key(k[0]), k[1], k[2]))
        rid = {k: i for i, k in enumerate(keys)}
        record_of = np.array([rid[(ch, int(s), int(e))] for ch, s, e in zip(chroms, lv["start"], lv["end"])],
                             dtype=np.uint32)
        recs = [record(ctx, tb, *k, fine=fine) for k in keys]
        # Décile SON parmi les intervalles du niveau : un score TSA-seq ne se lit qu'en relatif.
        sv = np.array([r["speckles"]["v"]["score"] if "v" in r["speckles"] else np.nan for r in recs], dtype=float)
        ok = np.isfinite(sv)
        if ok.any():
            edges = np.quantile(sv[ok], np.linspace(0, 1, 11))
            for r, x in zip(recs, sv):
                if math.isfinite(x):
                    r["speckles"]["v"]["decile"] = int(min(10, max(1, np.searchsorted(edges, x, side="right"))))
        ens = ens_n if name == "noyau" else ens_a if name == "apercu" else None
        model = model_fields(lv, ens, parent if fine else None, R)

        blocks = []
        for k0 in range(0, len(recs), BLOCK):
            part = recs[k0:k0 + BLOCK]
            spec = put(_blob(part))
            blocks.append({**spec, "first": k0, "n": len(part), "chrom": keys[k0][0]})
        ro = put(_deflate(record_of.tobytes()))
        mo = put(_blob(model))
        levels_meta.append({"name": name, "n_beads": int(len(record_of)), "records": len(recs), "fine": fine,
                            "record_of": {**ro, "dtype": "uint32"}, "model": mo, "blocks": blocks})
        raw = sum(len(json.dumps(r, ensure_ascii=False, separators=(",", ":"))) for r in recs)
        report["levels"][name] = {"records": len(recs), "beads": int(len(record_of)),
                                  "stored": sum(b["length"] for b in blocks) + ro["length"] + mo["length"],
                                  "raw_json": raw, "seconds": round(time.time() - t1, 1)}
        log(f"  {name:9s} {len(recs):5d} intervalles, {len(record_of):5d} billes, "
            f"{report['levels'][name]['stored'] / 1e6:.2f} Mo, {time.time() - t1:.0f} s")

    tables = {}
    for tname, parts in tb.blobs().items():
        tables[tname] = [{**meta, **put(_blob(content))} for meta, content in parts]
        report.setdefault("tables", {})[tname] = sum(t["length"] for t in tables[tname])
    log("  tables   " + " · ".join(f"{k} {v / 1e6:.2f} Mo" for k, v in report["tables"].items()))

    # Blocs de la chaîne hg38 → T2T sur les régions fines : le navigateur lève une base.
    chains = {"regions": [], "targets": []}
    seen: dict[int, int] = {}
    for chrom, s, e in fine_regions:
        bl = []
        for t0_, t1_, q0, k in ctx.hs1.blocks(chrom, s, e):
            if k not in seen:
                seen[k] = len(chains["targets"])
                chains["targets"].append([*ctx.hs1.targets[k], ctx.hs1.scores[k]])
            bl.append([t0_, t1_, q0, seen[k]])
        chains["regions"].append({"chrom": chrom, "start": s, "end": e, "blocks": bl})

    m = ctx.m
    used = sorted({"hg38.chrom.sizes", "chain.hg38.to.hs1", "hg38.cytoband", "rao2014.gm12878.subcompartments",
                   "rao2014.gm12878.domains", "rao2014.gm12878.loops", "lads.damid.lmnb1", "tsaseq.son",
                   "nads.hffc6.nucleolus", "nads.hffc6.control", *[k for k, _ in S.REPLI_PHASES],
                   "chain.hg19.to.hg38", "gencode.annotation.gtf", "gtex.median.tpm", "gnomad.constraint",
                   "encode.ccres", "chromhmm.gm12878", "ctcf.peaks.gm12878", "ctcf.motif.orientation",
                   "clinvar.vcf", "gwas.catalog", "mondo.obo", "hpo.obo", "hg38.2bit"})
    from ..export.build import pipeline_version

    header = {
        "format": "geno-fiche",
        "version": VERSION,
        "assembly": "GRCh38",
        "cell_type": "GM12878",
        "g3d": {"path": Path(g3d_path).name, "sha256": hashlib.sha256(Path(g3d_path).read_bytes()).hexdigest()},
        "pipeline": pipeline_version(),
        "evidence": EVIDENCE,
        "sections": SECTIONS,
        "fine_only": FINE_ONLY,
        "labels": LABELS,
        "sources": {k: {"url": m[k].url, "sha256": m[k].sha256, "license": m[k].license, "notes": m[k].notes}
                    for k in used},
        "tissues": ctx.tissues,
        "clinvar_classes": S.SIG,
        "repli_terciles": [_r(x, 1) for x in ctx.repli_terciles],
        "nuclear_radius_nm": R,
        "lifted": ctx.lift_report,
        "chains": {"hs1": chains},
        "levels": levels_meta,
        "tables": tables,
        "coordinates": "0-based, demi-ouvert ; affichage 1-based",
    }
    hbytes = _deflate(json.dumps(header, ensure_ascii=False, separators=(",", ":")).encode())
    pre = struct.pack("<8sII32s", MAGIC, VERSION, len(hbytes), hashlib.sha256(hbytes).digest()).ljust(PREAMBLE, b"\0")
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as f:
        f.write(pre)
        f.write(hbytes)
        for b in blobs:
            f.write(b)
    report.update({"path": str(out), "bytes": PREAMBLE + len(hbytes) + offset, "header": len(hbytes),
                   "seconds": round(time.time() - t0, 1)})
    return report


# ---------------------------------------------------------------------------
# Lecture (tests, `geno fiche check`, et pour garder le format honnête)
# ---------------------------------------------------------------------------

class Fiche:
    def __init__(self, path: str | Path):
        self.blob = Path(path).read_bytes()
        magic, version, hlen, hsha = struct.unpack_from("<8sII32s", self.blob)
        if magic != MAGIC:
            raise ValueError("ce n'est pas un fichier de fiches")
        if version != VERSION:
            raise ValueError(f"fiche version {version}, ce lecteur lit la {VERSION}")
        hb = self.blob[PREAMBLE:PREAMBLE + hlen]
        if hashlib.sha256(hb).digest() != hsha:
            raise ValueError("en-tête : empreinte divergente")
        self.header = json.loads(_inflate(hb))
        self.base = PREAMBLE + hlen
        self._cache: dict = {}

    def _get(self, spec: dict) -> bytes:
        d = self.blob[self.base + spec["offset"]: self.base + spec["offset"] + spec["length"]]
        if digest(d) != spec["sha256"]:
            raise ValueError("bloc : empreinte divergente")
        return _inflate(d)

    def level(self, name: str) -> dict:
        return next(lv for lv in self.header["levels"] if lv["name"] == name)

    def records(self, name: str) -> list[dict]:
        if ("r", name) not in self._cache:
            lv = self.level(name)
            self._cache[("r", name)] = [r for b in lv["blocks"] for r in json.loads(self._get(b))]
        return self._cache[("r", name)]

    def model(self, name: str) -> list[dict]:
        if ("m", name) not in self._cache:
            self._cache[("m", name)] = json.loads(self._get(self.level(name)["model"]))
        return self._cache[("m", name)]

    def record_of(self, name: str) -> np.ndarray:
        return np.frombuffer(self._get(self.level(name)["record_of"]), dtype=np.uint32)

    def table(self, name: str) -> list:
        if ("t", name) not in self._cache:
            rows: list = []
            for spec in self.header["tables"][name]:
                part = json.loads(self._get(spec))
                if name == "genes":
                    n = spec["n"]
                    part = [{"id": part["id"][i], "name": part["name"][i], "type": part["type"][i],
                             "strand": part["strand"][i], "start": part["start"][i], "end": part["end"][i],
                             "mane": part["mane"][i] == "1", "gtex": part["gtex"][i], "gnomad": part["gnomad"][i]}
                            for i in range(n)]
                rows += part
            self._cache[("t", name)] = rows
        return self._cache[("t", name)]

    def expand(self, rec: dict) -> dict:
        """Remet en place ce que les fiches partagent (gènes, maladies, traits)."""
        return expand(rec, self.table("genes"), self.table("diseases"), self.table("traits"))

    def card(self, name: str, bead_id: int) -> dict:
        """La fiche d'une bille (identifiant ≥ 1, comme le picking) : enregistrement + modèle."""
        k = bead_id - 1
        return {**self.expand(self.records(name)[int(self.record_of(name)[k])]), **self.model(name)[k]}


GENE_KEYS = ("id", "name", "type", "strand", "start", "end", "mane")


def expand(rec: dict, genes: list, diseases: list, traits: list) -> dict:
    """Un enregistrement aux indices remplacés par ce qu'ils désignent — la même opération
    que `FicheReader.card` côté navigateur."""
    r = dict(rec)
    g = dict(r["genes"])
    v = g["v"]
    g["v"] = {"genes": [{k: genes[i][k] for k in GENE_KEYS} for i in v["genes"]]}
    if "nearest" in v:
        g["v"]["nearest"] = [{**{k: genes[n["g"]][k] for k in GENE_KEYS}, "distance": n["distance"],
                              "side": n["side"]} for n in v["nearest"]]
    r["genes"] = g
    if "v" in r["expression"]:
        r["expression"] = {**r["expression"], "v": {
            genes[i]["id"]: dict(zip(("lcl", "tissue", "max"), genes[i]["gtex"])) for i in r["expression"]["v"]}}
    if "v" in r["contrainte"]:
        r["contrainte"] = {**r["contrainte"], "v": {
            genes[i]["id"]: dict(zip(("pli", "loeuf"), genes[i]["gnomad"])) for i in r["contrainte"]["v"]}}
    cv = dict(r["clinvar"]["v"])

    def var(x):
        return {**x, "diseases": [diseases[d] for d in x["diseases"]]}

    if "variants" in cv:
        cv["variants"] = [var(x) for x in cv["variants"]]
        cv["plp"] = sorted((x for x in cv["variants"] if x["sig"] == "plp"), key=lambda x: (-x["stars"], x["pos"]))
    else:
        cv["plp"] = [var(x) for x in cv["plp"]]
    r["clinvar"] = {**r["clinvar"], "v": cv}
    gw = dict(r["gwas"]["v"])
    gw["top"] = [{**{k: x for k, x in h.items() if k != "t"}, "trait": traits[h["t"]][0],
                  "mapped": traits[h["t"]][1], "uri": traits[h["t"]][2]} for h in gw["top"]]
    r["gwas"] = {**r["gwas"], "v": gw}
    return r


def missing(card: dict, header: dict, fine: bool) -> list[str]:
    """Les champs attendus qui n'ont ni valeur ni absence motivée. Vide = fiche complète."""
    out = []
    want = [f for fs in header["sections"].values() for f in fs]
    if fine:
        want += [f for fs in header["fine_only"].values() for f in fs]
    for f in want:
        d = card.get(f)
        if not isinstance(d, dict) or "e" not in d or "s" not in d:
            out.append(f)
        elif "v" not in d and not d.get("absent"):
            out.append(f)
        elif "v" in d and d["v"] is None:
            out.append(f)
    return out
