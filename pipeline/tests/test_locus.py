"""La fiche de locus : le critère de fin de la semaine 2.

`geno query chr7:5,527,000-5,530,600` doit rendre le locus ACTB complet.
"""

import pytest

from geno_pipeline.locus import RegionError, parse_region, render, report

ACTB = "chr7:5,527,000-5,530,600"


@pytest.mark.parametrize(
    "text,expected",
    [
        ("chr7:5,527,000-5,530,600", ("chr7", 5526999, 5530600)),
        ("chr7:5527000-5530600", ("chr7", 5526999, 5530600)),
        ("chr7:5_527_000-5_530_600", ("chr7", 5526999, 5530600)),
        ("  chr17:7668402-7687550  ", ("chr17", 7668401, 7687550)),
    ],
)
def test_region_parsing_is_one_based_inclusive(text, expected):
    assert parse_region(text) == expected


@pytest.mark.parametrize("bad", ["chr7", "chr7:100", "chr7:0-100", "chr7:500-100", "n'importe quoi"])
def test_bad_regions_are_refused(bad):
    with pytest.raises(RegionError):
        parse_region(bad)


def test_whole_chromosome_needs_sizes():
    assert parse_region("chr7", {"chr7": 159_345_973}) == ("chr7", 0, 159_345_973)


def test_actb_locus_is_complete(store):
    """Le critère de la feuille de route, vérifié piste par piste."""
    rep = report(store, ACTB)

    genes = rep.tracks["genes"]
    assert [g.attrs["name"] for g in genes] == ["ACTB"]
    assert genes[0].attrs["strand"] == "-"
    assert (genes[0].start, genes[0].end) == (5527150, 5530601)

    assert len(rep.tracks["exons"]) == 6
    assert [b.attrs["name"] for b in rep.tracks["cytoband"]] == ["7p22.1"]

    # ACTB est sur le brin moins : le promoteur est à la coordonnée la plus haute.
    classes = [c.attrs["name"] for c in rep.tracks["ccre"]]
    assert "PLS" in classes
    pls = next(c for c in rep.tracks["ccre"] if c.attrs["name"] == "PLS")
    assert pls.start > genes[0].start + (genes[0].end - genes[0].start) // 2

    states = [s.attrs["name"] for s in rep.tracks["chromhmm"]]
    assert "Tx" in states and "TssA" in states


def test_flank_reveals_the_convergent_ctcf_pair(store):
    """Sans marge, la piste CTCF d'ACTB est vide — et c'est exact.

    Les deux ancres d'une boucle convergente encadrent le gène : elles ne le
    chevauchent pas. Une fiche de locus qui ne regarde que la fenêtre stricte
    rate donc toujours la boucle du gène qu'elle décrit.
    """
    assert report(store, ACTB).tracks["ctcf"] == []

    rep = report(store, ACTB, flank=5_000)
    anchors = rep.tracks["ctcf"]
    assert [c.attrs["strand"] for c in anchors] == ["+", "-"]   # motifs convergents
    assert anchors[0].end < 5527150 and anchors[1].start > 5530601
    assert rep.flank == 5_000
    assert rep.span == 3601                                      # la région affichée ne bouge pas
    assert "± 5,000 pb" in render(rep)


def test_negative_flank_is_refused(store):
    with pytest.raises(RegionError):
        report(store, ACTB, flank=-1)


def test_report_carries_its_provenance(store):
    rep = report(store, ACTB)
    assert rep.source == "fixtures"          # une fiche de test se reconnaît
    assert rep.assembly == "GRCh38"
    assert "fixtures" in render(rep)


def test_track_filter(store):
    rep = report(store, ACTB, ["genes"])
    assert set(rep.tracks) == {"genes"}


def test_unknown_track_names_the_known_ones(store):
    with pytest.raises(KeyError, match="genes"):
        report(store, ACTB, ["inexistante"])


def test_empty_region_is_reported_as_empty(store):
    rep = report(store, "chr7:100,000,000-100,001,000")
    assert rep.is_empty()
    assert "aucune annotation" in render(rep)


def test_json_shape_is_stable(store):
    d = report(store, ACTB).to_dict()
    assert d["region"]["display"] == "chr7:5,527,000-5,530,600"
    assert d["region"]["span"] == 3601
    assert d["tracks"]["genes"][0]["name"] == "ACTB"
    assert d["tracks"]["genes"][0]["start"] == 5527150


def test_render_shows_every_track(store):
    text = render(report(store, ACTB))
    for label in ("bande cytogénétique", "gènes", "exons", "sites CTCF"):
        assert label in text
