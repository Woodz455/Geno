"""L'index recompté depuis les fichiers bruts, par un autre chemin — et les attendus de la
lecture croisée et de la mesure dans la page.

Autre chemin : une expression régulière sur les lignes brutes du VCF ClinVar, du TSV du
catalogue GWAS et du GTF de GENCODE, sans `fiche.sources` ni `_cached`. Le témoin décale chaque
position d'une base : il doit diverger partout.
"""

from __future__ import annotations

import io
import json
import random
import re
import zipfile
from pathlib import Path

from .build import MODEL, Index
from .resolve import resolve

_RS_INFO = re.compile(r"(?:^|;)RS=([\d|]+)")
_GENE = re.compile(r'gene_id "([^"]+)".*?gene_name "([^"]+)"')


def _raw_clinvar(path: Path, wanted: set[int]) -> dict[int, set]:
    import gzip

    out: dict[int, set] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            m = _RS_INFO.search(c[7])
            if not m:
                continue
            for r in m.group(1).split("|"):
                if r and int(r) in wanted:
                    chrom = c[0] if c[0].startswith("chr") else ("chrM" if c[0] == "MT" else "chr" + c[0])
                    out.setdefault(int(r), set()).add((chrom, int(c[1]) - 1))
    return out


def _raw_gwas(path: Path, wanted: set[int]) -> dict[int, set]:
    out: dict[int, set] = {}
    with zipfile.ZipFile(path) as z:
        name = next(n for n in z.namelist() if n.endswith(".tsv"))
        with z.open(name) as raw:
            f = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            col = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
            for line in f:
                c = line.rstrip("\n").split("\t")
                toks = re.findall(r"\brs(\d+)\b", c[col["SNPS"]])
                if len(toks) != 1 or int(toks[0]) not in wanted:
                    continue
                chrom, pos = c[col["CHR_ID"]], c[col["CHR_POS"]]
                if not chrom or not pos.isdigit():
                    continue
                try:
                    float(c[col["P-VALUE"]])
                except ValueError:
                    continue
                chrom = chrom if chrom.startswith("chr") else ("chrM" if chrom == "MT" else "chr" + chrom)
                out.setdefault(int(toks[0]), set()).add((chrom, int(pos) - 1))
    return out


def _raw_genes(path: Path, wanted: set[str]) -> dict[str, set]:
    import gzip

    out: dict[str, set] = {}
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    with (gzip.open(path, "rt") if gz else open(path)) as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 8)
            if c[2] != "gene":
                continue
            m = _GENE.search(c[8])
            if not m:
                continue
            gid, name = m.group(1).split(".")[0], m.group(2)
            for key in (name.upper(), gid.upper()):
                if key in wanted:
                    out.setdefault(key, set()).add((c[0], int(c[3]) - 1, int(c[4])))
    return out


def check(index_path: Path, m, n: int = 500, seed: int = 16, log=print) -> dict:
    from ..fiche import sources as S

    ix = Index(index_path)
    rng = random.Random(seed)
    rows = ix.all_rs()
    present = sorted(set(int(x) for x in rows["rs"]))
    sample = set(rng.sample(present, min(n, len(present)))) | {334}
    absent: set[int] = set()
    pres = set(present)
    while len(absent) < n:
        x = rng.randint(1, max(present))
        if x not in pres:
            absent.add(x)
    wanted = sample | absent
    log(f"  relecture brute : ClinVar, GWAS, pour {len(sample)} rsID présents et {len(absent)} absents …")
    raw = _raw_clinvar(S.locked(m, "clinvar.vcf").path, wanted)
    for k, v in _raw_gwas(S.locked(m, "gwas.catalog").path, wanted).items():
        raw.setdefault(k, set()).update(v)
    chroms = set(ix.chroms)
    eq = shifted = 0
    bad: list[int] = []
    for r in sample:
        got = {(c, p) for c, p, _ in ix.rs(r)}
        want = {x for x in raw.get(r, set()) if x[0] in chroms}
        eq += got == want
        shifted += {(c, p + 1) for c, p in got} == want
        if got != want:
            bad.append(r)
    absent_ok = sum(not ix.rs(r) and not raw.get(r) for r in absent)

    keys = _gene_keys(ix)
    gs = set(rng.sample(keys, min(n, len(keys)))) | {"BRCA1", "TP53", "HBB"}
    rawg = _raw_genes(S.locked(m, "gencode.annotation.gtf").path, gs)
    geq = sum({(r[3], r[4], r[5]) for r in ix.genes(k)} == rawg.get(k, set()) for k in gs)

    named = {
        "rs334": resolve(ix, "rs334"),
        "BRCA1": resolve(ix, "BRCA1"),
        "chr17:43,044,295": resolve(ix, "chr17:43,044,295"),
    }
    rep = {
        "rs_present": {"n": len(sample), "equal": eq, "shifted_control": shifted, "mismatch": sorted(bad)[:20]},
        "rs_absent": {"n": len(absent), "absent_both": absent_ok},
        "genes": {"n": len(gs), "equal": geq},
        "named": named,
    }
    log(f"  rsID présents : {eq}/{len(sample)} égaux aux fichiers bruts (témoin décalé : {shifted})")
    log(f"  rsID absents  : {absent_ok}/{len(absent)} absents des deux côtés")
    log(f"  gènes         : {geq}/{len(gs)} égaux au GTF brut")
    for q, r in named.items():
        h = r["hits"]
        log(f"  {q:18s} → " + "; ".join(f"{x['chrom']}:{x['start'] + 1:,}-{x['end']:,}" for x in h))
    return rep


def _gene_keys(ix: Index) -> list[str]:
    from ..export.g3d import _inflate

    keys: list[str] = []
    for b in ix.header["tables"]["genes"]["blocks"]:
        keys += [r[0] for r in json.loads(_inflate(ix._blob(b)))]
    return sorted(set(keys))


def dump(index_path: Path, out: Path, seed: int = 1616) -> dict:
    """Les requêtes des critères 1 et 2 de la semaine 16, avec la résolution Python : la page et
    le lecteur TypeScript doivent rendre la même chose."""
    from ..export.g3d import _inflate

    ix = Index(index_path)
    rng = random.Random(seed)
    model = list(MODEL)
    sizes = ix.sizes
    names = []
    for b in ix.header["tables"]["genes"]["blocks"]:
        names += [r for r in json.loads(_inflate(ix._blob(b))) if r[0] == r[1].upper() and r[3] in model
                  and r[7] == "protein_coding"]
    picked = rng.sample(names, 199)
    genes = ["BRCA1"] + [r[1] for r in picked]
    # Pour la grammaire : un nom sur dix en minuscules, un sur vingt par son identifiant ENSG.
    genes = [g.lower() if k % 10 == 9 else g for k, g in enumerate(genes)]
    for k in range(4, 200, 20):
        genes[k] = picked[k - 1][2]
    positions = ["chr17:43,044,295"]
    while len(positions) < 200:
        c = rng.choice(model)
        p = rng.randint(1, sizes[c])
        style = len(positions) % 4
        if style == 0:
            positions.append(f"{c}:{p:,}")
        elif style == 1:
            positions.append(f"{c[3:]}:{p}")
        elif style == 2:
            positions.append(f"{c} {p}")
        else:
            q = min(sizes[c], p + rng.randint(1, 200_000))
            positions.append(f"{c}:{p:,}-{q:,}")
    rows = ix.all_rs()
    model_idx = {ix.chroms.index(c) for c in model}
    inm = sorted({int(r) for r, c in zip(rows["rs"], rows["chrom"]) if int(c) in model_idx})
    rss = ["rs334"] + [f"rs{x}" for x in rng.sample(inm, 199)]
    pres = set(int(x) for x in rows["rs"])
    absent = []
    while len(absent) < 20:
        x = rng.randint(1, int(rows["rs"].max()))
        if x not in pres:
            absent.append(f"rs{x}")
    malformed = ["", "chr17:", "chr17:0", "chr17:20-10", "rs", "rs12x", "BRCA1 TP53", "%%%", "chr17:43,044,295:12",
                 "::", "-5", "12345", "chr1:abc", "rsrs1", "gène", "BRCA1;", "a" * 41, "chr1:1-", "x:y", "17-43"]
    other = ["chrY:2,786,855", "chrY:10,000,000", "chr23:5", "chr1:300,000,000", "chrM:100"]
    queries = [("gene", q) for q in genes] + [("position", q) for q in positions] + [("rs", q) for q in rss] + \
              [("rs-absent", q) for q in absent] + [("malformed", q) for q in malformed] + [("other", q) for q in other]
    out.parent.mkdir(parents=True, exist_ok=True)
    data = [{"group": g, "query": q, "expected": resolve(ix, q)} for g, q in queries]
    out.write_text(json.dumps({"queries": data, "seed": seed}, ensure_ascii=False))
    groups: dict[str, int] = {}
    for g, _ in queries:
        groups[g] = groups.get(g, 0) + 1
    return {"queries": len(data), "groups": groups}
