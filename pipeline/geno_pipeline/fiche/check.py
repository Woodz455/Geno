"""Ce que vaut le fichier de fiches, mesuré par un autre chemin que celui qui l'a écrit.

Trois mesures, chacune avec son témoin :

1. **Complétude** — toutes les billes de tous les niveaux. Un champ attendu sans valeur et
   sans absence motivée est un champ muet ; on en attend zéro. Témoin : la même mesure sur
   des fiches dont on a retiré un champ au hasard doit le trouver.
2. **Exactitude** — 200 billes tirées au hasard, recomposées sans `fiche.sources` :
   - gènes, cCRE, CTCF, cytobandes : le magasin d'intervalles de la semaine 2 (`intervals`,
     `parsers`), construit dans un dossier temporaire sur les fichiers verrouillés ;
   - ClinVar et GWAS : un balayage ligne à ligne des fichiers bruts, qui ne partage avec
     `sources` que le nom des colonnes.
   Témoin : les mêmes comparaisons, intervalles décalés d'une bille, doivent diverger.
3. **Faits connus** — écrits avant de lire le fichier.
"""

from __future__ import annotations

import gzip
import io
import random
import tempfile
import zipfile
from pathlib import Path

from ..intervals import Store, write_store
from ..parsers import parse_bed, parse_cytoband, parse_gtf
from . import sources as S
from .build import Fiche, missing


def completeness(f: Fiche) -> dict:
    out = {}
    for lv in f.header["levels"]:
        name = lv["name"]
        recs, model, ro = [f.expand(r) for r in f.records(name)], f.model(name), f.record_of(name)
        bad: dict[str, int] = {}
        absent: dict[str, int] = {}
        for k in range(lv["n_beads"]):
            card = {**recs[int(ro[k])], **model[k]}
            for fld in missing(card, f.header, lv["fine"]):
                bad[fld] = bad.get(fld, 0) + 1
            for fld, d in card.items():
                if "absent" in d:
                    absent[fld] = absent.get(fld, 0) + 1
        out[name] = {"beads": lv["n_beads"], "mute": bad, "absent": absent}
    # Témoin : un champ retiré d'une fiche doit être vu.
    lv = f.header["levels"][0]
    card = f.card(lv["name"], 1)
    card.pop("clinvar")
    out["control"] = missing(card, f.header, lv["fine"])
    return out


def _store(m: dict, root: Path) -> Store:
    def tracks():
        yield "genes", parse_gtf(S.locked(m, "gencode.annotation.gtf").path, "gene")
        yield "ccre", parse_bed(S.locked(m, "encode.ccres").path)
        yield "ctcf", parse_bed(S.locked(m, "ctcf.peaks.gm12878").path)
        yield "cytoband", (r for r in parse_cytoband(S.locked(m, "hg38.cytoband").path) if r[3]["band"])

    write_store(root, dict(tracks()), meta={"assembly": "GRCh38", "source": "manifeste verrouillé"})
    return Store(root)


def _scan_points(path: Path, intervals: list[tuple[str, int, int]], kind: str) -> list[dict]:
    """Compte, par intervalle, les lignes ClinVar (par classe) ou GWAS d'un fichier brut."""
    by: dict[str, list[tuple[int, int, int]]] = {}
    for k, (c, s, e) in enumerate(intervals):
        by.setdefault(c.removeprefix("chr"), []).append((s, e, k))
    out = [dict() for _ in intervals]

    def hit(chrom: str, pos0: int, label: str):
        for s, e, k in by.get(chrom, ()):
            if s <= pos0 < e:
                out[k][label] = out[k].get(label, 0) + 1

    if kind == "clinvar":
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                chrom, pos, _, _, _, _, _, info = line.split("\t", 8)[:8]
                if chrom not in by:
                    continue
                sig = ""
                for kv in info.split(";"):
                    if kv.startswith("CLNSIG="):
                        sig = kv[7:]
                hit(chrom, int(pos) - 1, S._sig(sig))
    else:
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith(".tsv"))
            with z.open(name) as raw:
                fh = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                head = fh.readline().rstrip("\n").split("\t")
                ci, pi, vi = head.index("CHR_ID"), head.index("CHR_POS"), head.index("P-VALUE")
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if c[ci] in by and c[pi].isdigit():
                        try:
                            float(c[vi])
                        except ValueError:
                            continue
                        hit(c[ci], int(c[pi]) - 1, "n")
    return out


def exactness(f: Fiche, n: int = 200, seed: int = 13, shift: bool = False) -> dict:
    m = S.manifest()
    rng = random.Random(seed)
    picks = []
    for lv in f.header["levels"]:
        recs = f.records(lv["name"])
        for k in rng.sample(range(len(recs)), min(len(recs), n // len(f.header["levels"]))):
            picks.append((lv["name"], k, f.expand(recs[k])))
    intervals = []
    for name, k, r in picks:
        h = r["hg38"]["v"]
        d = (h["end"] - h["start"]) if shift else 0
        intervals.append((h["chrom"], h["start"] + d, h["end"] + d))

    with tempfile.TemporaryDirectory() as tmp:
        store = _store(m, Path(tmp))
        res = {"beads": len(picks), "genes": 0, "ccre": 0, "ctcf": 0, "cytoband": 0, "clinvar": 0, "gwas": 0}
        cv = _scan_points(S.locked(m, "clinvar.vcf").path, intervals, "clinvar")
        gw = _scan_points(S.locked(m, "gwas.catalog").path, intervals, "gwas")
        first_diff = {}
        for (name, k, r), (chrom, s, e), cvc, gwc in zip(picks, intervals, cv, gw):
            q = dict(store.query(chrom, s, e, ["genes", "ccre", "ctcf", "cytoband"]))
            genes = sorted(g.attrs["gene_id"].split(".")[0] for g in q["genes"])
            mine = sorted(g["id"] for g in r["genes"]["v"]["genes"])
            ok = {
                "genes": genes == mine,
                "ccre": len(q["ccre"]) == r["ccre"]["v"]["n"] + r["ccre"]["v"]["inactive"],
                "ctcf": len(q["ctcf"]) == r["ctcf"]["v"]["n"],
                "cytoband": [b.attrs["name"] for b in q["cytoband"]] == (r["cytobande"].get("v") or []),
                "clinvar": {k_: v for k_, v in cvc.items() if v} ==
                           {k_: v for k_, v in r["clinvar"]["v"]["counts"].items() if v},
                "gwas": gwc.get("n", 0) == r["gwas"]["v"]["n"],
            }
            for key, good in ok.items():
                res[key] += int(good)
                if not good and key not in first_diff:
                    first_diff[key] = f"{name} {chrom}:{s}-{e}"
        res["first_diff"] = first_diff
    return res


def facts(f: Fiche) -> list[tuple[str, bool, str]]:
    """Ce qu'on sait sans le fichier, vérifié dans le fichier."""
    out = []
    fin = next(lv["name"] for lv in f.header["levels"] if lv["fine"] and "tp53" in lv["name"])
    recs = f.records(fin)
    r = f.expand(next(x for x in recs if x["hg38"]["v"]["start"] <= 7_673_801 < x["hg38"]["v"]["end"]))
    out.append(("TP53 c.818 en 17p13.1", r["cytobande"].get("v") == ["17p13.1"], str(r["cytobande"].get("v"))))
    g = {x["name"]: x for x in r["genes"]["v"]["genes"]}
    out.append(("le monomère de c.818 porte TP53, brin −, MANE", "TP53" in g and g["TP53"]["strand"] == "-"
                and g["TP53"]["mane"], str(g.get("TP53"))))
    c = r["contrainte"].get("v", {}).get("ENSG00000141510")
    out.append(("TP53 : pLI et LOEUF de gnomAD présents", c is not None, str(c)))
    v = [x for x in r["clinvar"]["v"]["variants"] if x["pos"] == 7_673_801]
    # R273H = c.818G>A sur le brin codant (−), donc C>T sur le brin + de hg38.
    rs = [x for x in v if x.get("rs") == "rs28934576" and x["alt"] == "T" and x["sig"] == "plp"]
    out.append(("rs28934576 C>T (R273H) pathogène à chr17:7 673 802", bool(rs),
                f"{len(v)} variants à cette base ; " + "; ".join(f"{x['id']} {x['clnsig']}" for x in v[:4])))
    t = r["t2t"].get("v")
    out.append(("T2T levé, même chromosome", bool(t) and t["chrom"] == "chr17", str(t)))
    fa = next(lv["name"] for lv in f.header["levels"] if lv["fine"] and "tp53" not in lv["name"])
    ra = f.expand(next(x for x in f.records(fa) if x["hg38"]["v"]["start"] <= 5_529_000 < x["hg38"]["v"]["end"]))
    out.append(("ACTB en 7p22.1", ra["cytobande"].get("v") == ["7p22.1"], str(ra["cytobande"].get("v"))))
    out.append(("ACTB : exprimé en LCL (> 100 TPM)",
                (ra["expression"].get("v", {}).get("ENSG00000075624", {}).get("lcl") or 0) > 100,
                str(ra["expression"].get("v", {}).get("ENSG00000075624"))))
    return out


def dump_expected(f: Fiche, out: Path, n: int = 40, seed: int = 7) -> None:
    """Ce que le lecteur Python lit, pour la lecture croisée TypeScript (`tools/test-fiche.mjs`) :
    des fiches tirées au hasard dans chaque niveau, et des positions levées en T2T."""
    import json

    m = S.manifest()
    rng = random.Random(seed)
    cards = []
    for lv in f.header["levels"]:
        for bead in sorted(rng.sample(range(1, lv["n_beads"] + 1), min(n, lv["n_beads"]))):
            cards.append({"level": lv["name"], "bead": bead, "card": f.card(lv["name"], bead)})
    hs1 = S.chain_hs1(m)
    lifts = []
    for reg in f.header["chains"]["hs1"]["regions"]:
        for _ in range(500):
            p = rng.randrange(reg["start"], reg["end"])
            hit = hs1.point(reg["chrom"], p)
            lifts.append({"chrom": reg["chrom"], "pos": p,
                          "hit": None if hit is None else {"chrom": hit.chrom, "pos": hit.pos, "strand": hit.strand}})
    Path(out).write_text(json.dumps({"cards": cards, "lifts": lifts}, ensure_ascii=False))


def liftover_kmers(n: int = 10_000, seed: int = 13, shift: int = 0) -> dict:
    """La chaîne hg38 → T2T contrôlée sans UCSC : pour `n` positions hg38 tirées au hasard
    (hors N), le 21-mer centré doit être le même dans hs1.2bit à la position levée. Les deux
    assemblages diffèrent d'environ une SNV par kb : on attend ~98 %, pas 100 %. Témoin :
    `shift` décale la position levée ; à 1 kb, on attend ~0 %."""
    import re

    from ..sequence import TwoBit, revcomp

    m = S.manifest()
    hs1 = S.chain_hs1(m)
    g = TwoBit(S.locked(m, "hg38.2bit").path)
    h = TwoBit(S.locked(m, "hs1.2bit").path)
    sizes = {c: g.length(c) for c in S.MAIN if c in g.offsets}
    names = [c for c in sizes if re.match(r"^chr([0-9]+|X|Y)$", c)]
    w = [sizes[c] for c in names]
    rng = random.Random(seed)
    tried = lifted = same = base = 0
    while tried < n:
        c = rng.choices(names, weights=w)[0]
        p = rng.randrange(20, sizes[c] - 20)
        s = g.fetch(c, p - 10, p + 11)
        if "N" in s:
            continue
        tried += 1
        hit = hs1.point(c, p)
        if hit is None:
            continue
        lifted += 1
        q = hit.pos + shift
        if q - 10 < 0 or q + 11 > h.length(hit.chrom):
            continue
        t = h.fetch(hit.chrom, q - 10, q + 11)
        if hit.strand == "-":
            t = revcomp(t)
        same += t == s
        base += t[10] == s[10]
    return {"tried": tried, "lifted": lifted, "kmer": same, "base": base}
