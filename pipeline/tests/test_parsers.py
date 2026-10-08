"""La conversion de coordonnées est le seul endroit où une erreur d'une base passe inaperçue.

GTF est 1-based bornes incluses, BED est 0-based demi-ouvert. Tout le reste du
projet suppose la seconde convention. Ces tests verrouillent la frontière.
"""

from geno_pipeline.parsers import parse_bed, parse_cytoband, parse_gtf


def test_gtf_converts_to_zero_based_half_open(fixtures_dir):
    genes = {g[3]["name"]: g for g in parse_gtf(fixtures_dir / "gencode.fixture.gtf", "gene")}
    chrom, start, end, attrs = genes["ACTB"]
    # Le GTF porte 5527151-5530601 en 1-based inclus.
    assert (chrom, start, end) == ("chr7", 5527150, 5530601)
    assert end - start == 3451
    assert attrs["strand"] == "-"
    assert attrs["biotype"] == "protein_coding"
    assert attrs["gene_id"].startswith("ENSG00000075624")


def test_gtf_feature_filter_is_exact(fixtures_dir):
    path = fixtures_dir / "gencode.fixture.gtf"
    genes = list(parse_gtf(path, "gene"))
    exons = list(parse_gtf(path, "exon"))
    assert {g[3]["name"] for g in genes} == {"ACTB", "HBB", "TP53", "BRCA1"}
    assert sum(1 for e in exons if e[3]["name"] == "ACTB") == 6
    assert sum(1 for e in exons if e[3]["name"] == "HBB") == 3


def test_exons_stay_inside_their_gene(fixtures_dir):
    path = fixtures_dir / "gencode.fixture.gtf"
    bounds = {g[3]["name"]: (g[1], g[2]) for g in parse_gtf(path, "gene")}
    for chrom, start, end, attrs in parse_gtf(path, "exon"):
        lo, hi = bounds[attrs["name"]]
        assert lo <= start < end <= hi, f"exon hors du gène {attrs['name']}"


def test_bed_is_left_untouched(fixtures_dir):
    rows = list(parse_bed(fixtures_dir / "ccre.fixture.bed"))
    assert rows[0][:3] == ("chr7", 5530450, 5530800)
    assert rows[0][3]["name"] == "PLS"
    assert rows[0][3]["score"] == 956


def test_bed_strand_carries_motif_orientation(fixtures_dir):
    """Sans l'orientation du motif CTCF, la règle des motifs convergents est inapplicable."""
    rows = [r for r in parse_bed(fixtures_dir / "ctcf.fixture.bed") if r[0] == "chr7"]
    strands = [r[3]["strand"] for r in rows]
    assert strands[:2] == ["+", "-"]           # paire convergente encadrant ACTB
    assert all("strand" in r[3] for r in rows)


def test_cytoband_names_carry_the_chromosome(fixtures_dir):
    bands = {b[3]["name"]: b for b in parse_cytoband(fixtures_dir / "cytoband.fixture.txt")}
    chrom, start, end, attrs = bands["7p22.1"]
    assert (chrom, start, end) == ("chr7", 5000000, 7200000)
    assert attrs["arm"] == "p"
    assert attrs["stain"] == "gneg"
    assert "17q21.31" in bands


def test_comment_and_header_lines_are_skipped(fixtures_dir):
    # Les fixtures commencent toutes par des lignes # ; aucune ne doit être parsée.
    for row in parse_bed(fixtures_dir / "chromhmm.fixture.bed"):
        assert not row[0].startswith("#")
