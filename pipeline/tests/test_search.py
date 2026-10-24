"""Recherche (semaine 16) : la grammaire, et l'index écrit puis relu, contre des oracles écrits
à la main. Le vrai index est recompté depuis les fichiers bruts par `geno search check`."""

from __future__ import annotations

import pytest

from geno_pipeline.search import build as B
from geno_pipeline.search.query import parse
from geno_pipeline.search.resolve import resolve


@pytest.mark.parametrize("text,kind,chrom,start,end", [
    ("chr17:43,044,295", "position", "chr17", 43_044_294, 43_044_295),
    ("17:43044295", "position", "chr17", 43_044_294, 43_044_295),
    ("chrx 1234", "position", "chrX", 1233, 1234),
    ("MT:16", "position", "chrM", 15, 16),
    ("chr7:5,527,000-5,530,600", "interval", "chr7", 5_526_999, 5_530_600),
    ("chr7:5 527 000 – 5 530 600", "interval", "chr7", 5_526_999, 5_530_600),
])
def test_positions_are_one_based_inclusive_in_and_half_open_out(text, kind, chrom, start, end):
    q = parse(text)
    assert (q.kind, q.chrom, q.start, q.end) == (kind, chrom, start, end)


def test_rs_genes_ids_and_refusals():
    assert parse("RS334").kind == "rs" and parse("rs0334").key == "334"
    assert parse("brca1").kind == "gene" and parse("brca1").key == "BRCA1"
    assert parse("ENSG00000012048.25").key == "ENSG00000012048"
    for bad in ("", "chr17:0", "chr17:20-10", "12345", "BRCA1 TP53", "gène", "a" * 41):
        assert parse(bad).kind == "malformed", bad


def _index(tmp_path, monkeypatch, genes, rs_rows):
    """Un index écrit par `build`, sur des sources remplacées par des listes."""
    monkeypatch.setattr(B, "gene_rows", lambda m: sorted(
        [[n.upper(), n, i, c, s, e, "+", "protein_coding"] for n, i, c, s, e in genes]
        + [[i.upper(), n, i, c, s, e, "+", "protein_coding"] for n, i, c, s, e in genes],
        key=lambda r: (r[0], r[3], r[4], r[2])))
    monkeypatch.setattr(B, "clinvar_rs", lambda m: [(a, b, c) for a, b, c, s in rs_rows if s & 1])
    monkeypatch.setattr(B, "gwas_rs", lambda m: [(a, b, c) for a, b, c, s in rs_rows if s & 2])

    class E:
        def __init__(self, p):
            self.path, self.sha256, self.url = p, "0" * 64, "test"

    sizes = tmp_path / "sizes"
    sizes.write_text("chr1\t1000000\nchr11\t2000000\nchrY\t500000\nchr1_alt\t10\n")
    m = {k: E(sizes) for k in ("hg38.chrom.sizes", "gencode.annotation.gtf", "clinvar.vcf", "gwas.catalog")}
    from geno_pipeline.fiche import sources as S
    monkeypatch.setattr(S, "locked", lambda m, k: m[k])
    out = tmp_path / "t.search"
    B.build(m, out, log=lambda *a: None)
    return B.Index(out)


def test_a_name_spanning_two_blocks_returns_all_its_rows(tmp_path, monkeypatch):
    monkeypatch.setattr(B, "GENE_BLOCK", 4)
    genes = [("AAA", f"ENSG{k:011d}", "chr1", 10 * k, 10 * k + 5) for k in range(7)]     # 7 homonymes
    genes += [("ZZZ", "ENSG99999999999", "chr11", 5, 9)]
    ix = _index(tmp_path, monkeypatch, genes, [(1, "chr1", 0, 1)])
    assert len(ix.header["tables"]["genes"]["blocks"]) >= 4
    assert [r[4] for r in ix.genes("AAA")] == [0, 10, 20, 30, 40, 50, 60]
    assert ix.genes("AAB") == [] and len(ix.genes("ZZZ")) == 1


def test_rs_sources_merge_blocks_split_and_64_bit_numbers(tmp_path, monkeypatch):
    monkeypatch.setattr(B, "RS_BLOCK", 3)
    rows = [(5, "chr1", 100, 1), (5, "chr1", 100, 2), (5, "chr11", 7, 2), (6, "chr1", 1, 1),
            (7, "chr1", 2, 1), (8, "chr1", 3, 1), (9, "chr1", 4, 2), (8_766_588_714, "chr1", 9, 2),
            (11, "chr1_alt", 1, 1)]
    ix = _index(tmp_path, monkeypatch, [("A", "ENSG00000000001", "chr1", 0, 1)], rows)
    assert ix.rs(5) == [("chr1", 100, 3), ("chr11", 7, 2)]        # même position : sources fusionnées
    assert ix.rs(8_766_588_714) == [("chr1", 9, 2)]                 # au-delà de 2³² : bloc à part
    assert ix.rs(10) == [] and ix.rs(11) == []                      # contig alternatif : hors index
    blocks = ix.header["tables"]["rs"]["blocks"]
    assert all(b["last"] < c["first"] for b, c in zip(blocks, blocks[1:]))


def test_resolve_flags_out_of_model_and_out_of_range(tmp_path, monkeypatch):
    ix = _index(tmp_path, monkeypatch, [("SRY", "ENSG00000184895", "chrY", 100, 200)], [(334, "chr11", 1000, 1)])
    r = resolve(ix, "SRY")
    assert r["error"] == "" and r["hits"][0]["in_model"] is False
    assert resolve(ix, "rs334")["hits"][0] == {"label": "rs334", "chrom": "chr11", "start": 1000, "end": 1001,
                                               "sources": ["clinvar"], "in_model": True}
    assert resolve(ix, "chr1:1,000,001")["error"] == "out-of-range"
    assert resolve(ix, "chr22:5")["error"] == "unknown-chrom"
    assert resolve(ix, "rs335")["error"] == "unknown-rs"
    # L'en-tête est signé : un octet de plus et la lecture refuse.
    blob = bytearray(ix.blob)
    blob[B.PREAMBLE + 3] ^= 1
    p = tmp_path / "bad.search"
    p.write_bytes(bytes(blob))
    with pytest.raises(ValueError):
        B.Index(p)
