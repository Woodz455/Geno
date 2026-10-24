"""Une requête → ses cibles : la même réponse que `SearchIndex.resolve` en TypeScript."""

from __future__ import annotations

from .build import CLINVAR, GWAS, Index
from .query import parse


def resolve(ix: Index, text: str) -> dict:
    """{kind, hits: [{label, chrom, start, end, in_model, …}], error}. `error` vide si au moins
    une cible ; une cible hors du modèle (chrY, chrM : le modèle est 46,XX) reste une cible,
    marquée `in_model: false` — la page dit pourquoi elle ne vole pas."""
    q = parse(text)
    model = set(ix.header["model_chroms"])
    hits: list[dict] = []
    if q.kind == "malformed":
        return {"kind": q.kind, "hits": [], "error": "malformed"}
    if q.kind in ("gene", "gene-id"):
        for r in ix.genes(q.key):
            _, name, gid, chrom, s, e, strand, typ = r
            hits.append({"label": name, "id": gid, "chrom": chrom, "start": s, "end": e,
                         "strand": strand, "type": typ})
        err = "" if hits else "unknown-gene"
    elif q.kind == "rs":
        for chrom, pos, src in ix.rs(int(q.key)):
            hits.append({"label": f"rs{q.key}", "chrom": chrom, "start": pos, "end": pos + 1,
                         "sources": [n for n, b in (("clinvar", CLINVAR), ("gwas", GWAS)) if src & b]})
        err = "" if hits else "unknown-rs"
    else:
        size = ix.sizes.get(q.chrom)
        if size is None:
            return {"kind": q.kind, "hits": [], "error": "unknown-chrom"}
        if q.end > size:
            return {"kind": q.kind, "hits": [], "error": "out-of-range"}
        label = f"{q.chrom}:{q.start + 1:,}" if q.kind == "position" else f"{q.chrom}:{q.start + 1:,}-{q.end:,}"
        hits.append({"label": label, "chrom": q.chrom, "start": q.start, "end": q.end})
        err = ""
    for h in hits:
        h["in_model"] = h["chrom"] in model
    return {"kind": q.kind, "hits": hits, "error": err}
