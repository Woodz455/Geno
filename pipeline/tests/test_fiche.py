"""La fiche (semaine 13) : ce qui se teste sans les sources réelles.

Le fichier lui-même est vérifié sur les vraies données par `geno fiche check` (complétude,
exactitude recomposée par un autre chemin, faits connus) : ici, les briques, chacune contre
un oracle écrit à la main ou une force brute.
"""

from __future__ import annotations

import json
import random

import numpy as np
import pytest

from geno_pipeline.fiche import sources as S
from geno_pipeline.fiche.build import F, expand, missing

# MA0139.1, comptes JASPAR (19 positions, A C G T) — recopiés pour ne pas lire le réseau.
PFM = {
    "A": [87, 167, 281, 56, 8, 744, 40, 107, 851, 5, 333, 54, 12, 56, 104, 372, 82, 117, 402],
    "C": [291, 145, 49, 800, 903, 13, 528, 433, 11, 0, 3, 12, 0, 8, 733, 13, 482, 322, 181],
    "G": [76, 414, 449, 21, 0, 65, 334, 48, 32, 903, 566, 504, 890, 775, 5, 507, 307, 73, 266],
    "T": [459, 187, 134, 36, 2, 91, 11, 324, 18, 3, 9, 341, 8, 71, 67, 17, 37, 396, 59],
}


def _pwm():
    counts = np.array([PFM[b] for b in "ACGT"], dtype=float).T + 0.8
    prob = counts / counts.sum(axis=1, keepdims=True)
    return np.log2(prob / 0.25)


def _consensus():
    return "".join("ACGT"[int(np.argmax([PFM[b][i] for b in "ACGT"]))] for i in range(19))


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_motif_on_the_plus_strand_is_found_with_its_offset():
    rng = random.Random(1)
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    seq = flank(40) + _consensus() + flank(60)
    score, off, strand = S.motif_scan(seq, _pwm())
    assert (strand, off) == ("+", 40) and score > 0.95


def test_motif_on_the_minus_strand_is_oriented_minus():
    rng = random.Random(2)
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    seq = flank(25) + _revcomp(_consensus()) + flank(70)
    score, off, strand = S.motif_scan(seq, _pwm())
    assert (strand, off) == ("-", 25) and score > 0.95


def test_random_sequence_rarely_passes_the_motif_threshold():
    rng = random.Random(3)
    hits = sum(S.motif_scan("".join(rng.choice("ACGT") for _ in range(200)), _pwm())[0] >= S.MOTIF_MIN
               for _ in range(200))
    assert hits < 40                     # témoin : le seuil de 80 % n'est pas « tout »


def test_n_never_scores_as_a_motif():
    assert S.motif_scan("N" * 50, _pwm())[0] < 0.5


def test_intervals_query_matches_brute_force():
    rng = random.Random(4)
    rows = []
    for k in range(400):
        s = rng.randrange(0, 100_000)
        rows.append(("chr1", s, s + rng.randrange(1, 5_000), k))
    iv = S.Intervals.build(rows)
    for _ in range(300):
        a = rng.randrange(0, 100_000)
        b = a + rng.randrange(1, 8_000)
        got = sorted(iv.get("chr1", a, b))
        want = sorted(k for _, s, e, k in rows if s < b and e > a)
        assert got == want
    assert iv.get("chr2", 0, 10) == []


def test_overlap_bp_is_clipped_to_the_window():
    iv = S.Intervals.build([("chr1", 100, 200, "x"), ("chr1", 150, 400, "y")])
    assert sorted(iv.overlap_bp("chr1", 180, 300)) == [("x", 20), ("y", 120)]


def test_points_range_is_half_open():
    pts = S.Points.build([("chr1", p, p) for p in (5, 10, 10, 20)])
    assert sorted(pts.items[i] for i in pts.range("chr1", 10, 20)) == [10, 10]
    assert list(pts.range("chr1", 21, 30)) == []


@pytest.mark.parametrize("raw,code", [
    ("Pathogenic", "plp"), ("Likely_pathogenic", "plp"), ("Pathogenic/Likely_pathogenic", "plp"),
    ("Uncertain_significance", "vus"), ("Benign", "blb"), ("Likely_benign", "blb"),
    ("Conflicting_classifications_of_pathogenicity", "conf"), ("risk_factor", "other"),
    ("drug_response", "other"), ("", "other"),
])
def test_clinvar_classes(raw, code):
    assert S._sig(raw) == code


def _header():
    return {"sections": {"a": ["x", "y"]}, "fine_only": {"a": ["z"]}}


def test_a_field_with_a_value_or_a_reason_is_complete():
    card = {"x": F("mesure", "s", 1), "y": F("autre", "s", absent="pas mesuré ici")}
    assert missing(card, _header(), fine=False) == []


def test_a_mute_field_is_reported():
    card = {"x": F("mesure", "s", 1), "y": {"e": "autre", "s": "s"}}
    assert missing(card, _header(), fine=False) == ["y"]
    assert missing({"x": F("mesure", "s", None), "y": F("autre", "s", 2)}, _header(), fine=False) == ["x"]


def test_fine_levels_expect_their_extra_fields():
    card = {"x": F("mesure", "s", 1), "y": F("autre", "s", 2)}
    assert missing(card, _header(), fine=True) == ["z"]


def _rec(fine: bool):
    v = {"counts": {"plp": 2}, "plp_total": 2}
    variants = [{"sig": "plp", "stars": 1, "pos": 30, "diseases": [1]},
                {"sig": "vus", "stars": 2, "pos": 10, "diseases": []},
                {"sig": "plp", "stars": 3, "pos": 20, "diseases": [0, 1]}]
    if fine:
        v["variants"] = variants
    else:
        v["plp"] = [variants[2], variants[0]]
    return {
        "genes": F("annotation", "g", {"genes": [1], "nearest": [{"g": 0, "distance": 5, "side": "amont"}]}),
        "expression": F("autre", "x", [1]),
        "contrainte": F("annotation", "c", absent="non codant"),
        "clinvar": F("annotation", "cv", v),
        "gwas": F("annotation", "gw", {"n": 1, "top": [{"rs": "rs1", "pos": 3, "t": 0, "p": 1e-9, "study": "GCST1"}]}),
    }


GENES = [
    {"id": "G0", "name": "A", "type": "protein_coding", "strand": "+", "start": 0, "end": 5, "mane": True,
     "gtex": [1.0, 3, 9.0], "gnomad": [0.1, 1.2]},
    {"id": "G1", "name": "B", "type": "lncRNA", "strand": "-", "start": 9, "end": 50, "mane": False,
     "gtex": [2.5, 0, 2.5], "gnomad": None},
]
DIS = [[["MONDO:1"], "maladie un"], [["OMIM:2"], "maladie deux"]]
TRAITS = [["taille", "body height", "http://www.ebi.ac.uk/efo/EFO_0004339"]]


def test_expand_puts_the_shared_tables_back():
    r = expand(_rec(False), GENES, DIS, TRAITS)
    assert [g["name"] for g in r["genes"]["v"]["genes"]] == ["B"]
    assert r["genes"]["v"]["nearest"][0] | {} == {**{k: GENES[0][k] for k in ("id", "name", "type", "strand", "start",
                                                                               "end", "mane")},
                                                   "distance": 5, "side": "amont"}
    assert r["expression"]["v"] == {"G1": {"lcl": 2.5, "tissue": 0, "max": 2.5}}
    assert r["contrainte"]["absent"] == "non codant"
    assert r["clinvar"]["v"]["plp"][0]["diseases"] == [DIS[0], DIS[1]]
    assert r["gwas"]["v"]["top"][0]["trait"] == "taille" and "t" not in r["gwas"]["v"]["top"][0]


def test_fine_records_derive_their_plp_list_sorted_by_review_then_position():
    r = expand(_rec(True), GENES, DIS, TRAITS)
    assert [x["pos"] for x in r["clinvar"]["v"]["plp"]] == [20, 30]
    assert len(r["clinvar"]["v"]["variants"]) == 3


def test_expand_is_what_the_browser_reads_json_for_json():
    # Le navigateur relit du JSON : l'expansion ne doit rien produire qui ne s'y sérialise pas.
    json.dumps(expand(_rec(True), GENES, DIS, TRAITS), ensure_ascii=False)
