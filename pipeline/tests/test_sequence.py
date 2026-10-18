"""Semaine 12 — la séquence de référence et la base visée.

Le lecteur `.2bit` est vérifié sur des fichiers écrits ici, sans réseau. Le GTF aussi :
un transcrit sur le brin moins, des CDS dans plusieurs exons, et la règle du sens de
lecture — c'est exactement là qu'une coordonnée de codon se trompe.
"""

from __future__ import annotations

import gzip

import numpy as np
import pytest

from geno_pipeline.sequence import TwoBit, read_transcript, revcomp, write_twobit


def test_twobit_reads_back_any_interval_including_gaps(tmp_path):
    rng = np.random.default_rng(3)
    a = "".join(np.array(list("ACGT"))[rng.integers(0, 4, 10_007)])
    a = a[:500] + "N" * 37 + a[537:]
    b = "GATTACA" * 11
    write_twobit(tmp_path / "t.2bit", {"chrA": a, "chrB": b})
    with TwoBit(tmp_path / "t.2bit") as tb:
        assert tb.length("chrA") == len(a)
        for s, e in [(0, 1), (3, 9), (490, 560), (10_000, 10_007), (0, 10_007)]:
            assert tb.fetch("chrA", s, e) == a[s:e], (s, e)
        assert tb.fetch("chrB", 5, 40) == b[5:40]
        with pytest.raises(ValueError):
            tb.fetch("chrB", 0, len(b) + 1)
        with pytest.raises(KeyError):
            tb.fetch("chrZ", 0, 1)


def test_a_file_that_is_not_twobit_is_refused(tmp_path):
    (tmp_path / "x.2bit").write_bytes(b"\0" * 32)
    with pytest.raises(ValueError):
        TwoBit(tmp_path / "x.2bit")


GTF = """\
chrT\tTEST\tgene\t100\t400\t.\t-\t.\tgene_id "G1"; gene_name "TOY";
chrT\tTEST\ttranscript\t100\t400\t.\t-\t.\tgene_id "G1"; transcript_id "ENST1.4"; gene_name "TOY"; tag "MANE_Select";
chrT\tTEST\texon\t100\t150\t.\t-\t.\tgene_id "G1"; transcript_id "ENST1.4"; gene_name "TOY";
chrT\tTEST\tCDS\t120\t150\t.\t-\t0\tgene_id "G1"; transcript_id "ENST1.4"; gene_name "TOY";
chrT\tTEST\texon\t300\t400\t.\t-\t.\tgene_id "G1"; transcript_id "ENST1.4"; gene_name "TOY";
chrT\tTEST\tCDS\t300\t350\t.\t-\t0\tgene_id "G1"; transcript_id "ENST1.4"; gene_name "TOY";
chrT\tTEST\ttranscript\t100\t400\t.\t-\t.\tgene_id "G1"; transcript_id "ENST10.1"; gene_name "TOY";
"""


def test_coding_positions_run_against_the_genome_on_the_minus_strand(tmp_path):
    p = tmp_path / "t.gtf.gz"
    with gzip.open(p, "wt") as f:
        f.write(GTF)
    t = read_transcript(p, "ENST1")
    assert (t.gene, t.strand, t.transcript) == ("TOY", "-", "ENST1.4")
    assert "MANE_Select" in t.tags
    assert t.cds_length == 51 + 31
    # c.1 est l'extrémité 3' génomique du premier bloc codant dans le sens du transcrit.
    assert t.coding_position(1) == 350
    assert t.coding_position(51) == 300
    assert t.coding_position(52) == 150, "le saut d'intron se fait dans le sens du transcrit"
    assert t.codon(17) == (302, 301, 300)
    assert t.codon(18) == (150, 149, 148)
    with pytest.raises(ValueError):
        t.coding_position(83)


def test_a_transcript_id_prefix_is_not_a_match(tmp_path):
    p = tmp_path / "t.gtf"
    p.write_text(GTF)
    assert read_transcript(p, "ENST10").transcript == "ENST10.1"
    with pytest.raises(KeyError):
        read_transcript(p, "ENST2")


def test_revcomp():
    assert revcomp("ACGTN") == "NACGT"
