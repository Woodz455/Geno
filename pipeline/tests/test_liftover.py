"""Liftover par chaînes : un fichier écrit à la main, dont on connaît chaque réponse."""

from __future__ import annotations

import gzip

import numpy as np

from geno_pipeline.liftover import Chain, ucsc_name

# Chaîne 1 (score 1000) : chrA[100, 400) → chrB, brin +, avec un trou de 50 dans la source
# et de 20 dans la cible :  blocs [100,200)→[1000,1100), [250,400)→[1120,1270).
# Chaîne 2 (score 10)   : chrA[150, 180) → chrC (taille 500), brin − : bloc [150,180)→[10,40)
#                         compté depuis la fin de chrC, donc [460, 490) sur le brin +.
# Chaîne 3 (score 500)  : chrA[600, 700) → chrB, brin +, un seul bloc → [2000, 2100).
CHAINS = """chain 1000 chrA 1000 + 100 400 chrB 5000 + 1000 1270 1
100 50 20
150

chain 10 chrA 1000 + 150 180 chrC 500 - 10 40 2
30

chain 500 1 1000 + 600 700 chrB 5000 + 2000 2100 3
100
"""


def _chain(tmp_path, gz=False):
    p = tmp_path / "toy.chain"
    if gz:
        with gzip.open(p, "wt") as f:
            f.write(CHAINS)
    else:
        p.write_text(CHAINS)
    return Chain(p)


def test_points_inside_blocks_follow_their_offset(tmp_path):
    c = _chain(tmp_path)
    assert c.point("chrA", 100).pos == 1000
    assert c.point("chrA", 199).pos == 1099
    assert c.point("chrA", 250).pos == 1120
    assert c.point("chrA", 399).pos == 1269
    assert c.point("chrA", 100).chrom == "chrB"


def test_gaps_and_outside_do_not_lift(tmp_path):
    c = _chain(tmp_path)
    for p in (99, 200, 249, 400, 599, 700):
        assert c.point("chrA", p) is None, p
    assert c.point("chrZ", 150) is None


def test_overlapping_chains_keep_the_best_score(tmp_path):
    c = _chain(tmp_path)
    hit = c.point("chrA", 160)
    assert hit.chrom == "chrB" and hit.pos == 1060 and hit.strand == "+"


def test_minus_strand_counts_from_the_end_of_the_target(tmp_path):
    # La chaîne 2 seule, pour lire son brin − sans que la chaîne 1 gagne.
    p = tmp_path / "minus.chain"
    p.write_text("chain 10 chrA 1000 + 150 180 chrC 500 - 10 40 2\n30\n")
    c = Chain(p)
    hit = c.point("chrA", 150)
    assert hit.chrom == "chrC" and hit.strand == "-"
    assert hit.pos == 500 - 1 - 10
    assert c.point("chrA", 179).pos == 500 - 1 - 39


def test_ensembl_names_become_ucsc_names(tmp_path):
    c = _chain(tmp_path)
    assert "chr1" in c.chroms()
    assert c.point("chr1", 650).pos == 2050
    assert ucsc_name("MT") == "chrM" and ucsc_name("X") == "chrX" and ucsc_name("chr7") == "chr7"


def test_gzip_is_sniffed_not_guessed_from_the_suffix(tmp_path):
    c = _chain(tmp_path, gz=True)
    assert c.point("chrA", 100).pos == 1000


def test_interval_reports_coverage_and_skips_a_gap_at_its_ends(tmp_path):
    c = _chain(tmp_path)
    iv = c.interval("chrA", 100, 400)
    assert (iv.chrom, iv.start, iv.end, iv.split) == ("chrB", 1000, 1270, False)
    assert np.isclose(iv.mapped, 250 / 300)
    # Un bout dans le trou : on prend le premier point levable vers l'intérieur.
    iv = c.interval("chrA", 210, 300)
    assert (iv.start, iv.end) == (1120, 1170)


def test_interval_across_chains_is_flagged_split(tmp_path):
    p = tmp_path / "split.chain"
    p.write_text("chain 9 chrA 1000 + 0 100 chrB 5000 + 0 100 1\n100\n\n"
                 "chain 8 chrA 1000 + 100 200 chrD 900 + 0 100 2\n100\n")
    iv = Chain(p).interval("chrA", 50, 150)
    assert iv.split and iv.chrom == "chrB"
    assert not Chain(p).interval("chrA", 10, 90).split


def test_interval_with_no_block_does_not_lift(tmp_path):
    c = _chain(tmp_path)
    assert c.interval("chrA", 420, 580) is None


def test_blocks_for_embedding_are_the_ones_touching_the_window(tmp_path):
    c = _chain(tmp_path)
    b = c.blocks("chrA", 170, 260)
    assert [(t0, t1) for t0, t1, _, _ in b] == [(100, 200), (150, 180), (250, 400)]
    assert [(t0, t1) for t0, t1, _, _ in c.blocks("chrA", 190, 260)] == [(100, 200), (250, 400)]
