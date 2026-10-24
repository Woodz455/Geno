"""L'index doit donner exactement ce que donne une boucle naïve — jamais moins, jamais plus.

Un index d'intervalles qui perd silencieusement un chevauchement produit une fiche
de locus incomplète, ce qui est pire qu'une erreur : personne ne le remarque.
"""

import numpy as np
import pytest

from geno_pipeline.intervals import BLOCK, Store, Track, TrackWriter, _blockmax, _overlaps


def brute(starts, ends, qs, qe):
    return [i for i in range(len(starts)) if starts[i] < qe and ends[i] > qs]


def build(starts, ends):
    s = np.asarray(starts, dtype=np.int64)
    e = np.asarray(ends, dtype=np.int64)
    return s, e, _blockmax(e)


def test_half_open_boundaries():
    s, e, bm = build([100], [200])
    assert _overlaps(s, e, bm, 200, 300).tolist() == []   # touche à droite : pas de chevauchement
    assert _overlaps(s, e, bm, 0, 100).tolist() == []     # touche à gauche : pas de chevauchement
    assert _overlaps(s, e, bm, 199, 200).tolist() == [0]  # une base commune
    assert _overlaps(s, e, bm, 100, 101).tolist() == [0]


def test_degenerate_queries():
    s, e, bm = build([10, 20], [30, 40])
    assert _overlaps(s, e, bm, 25, 25).tolist() == []   # requête vide
    assert _overlaps(s, e, bm, 30, 20).tolist() == []   # requête inversée
    empty = np.empty(0, dtype=np.int64)
    assert _overlaps(empty, empty, _blockmax(empty), 0, 100).tolist() == []


def test_long_interval_is_not_missed():
    """Le cas qui casse un index naïf : un gène très long tout au début.

    CNTNAP2 fait 2,3 Mb. Si l'index ne remonte pas assez loin, il disparaît de
    toutes les requêtes situées après lui.
    """
    starts = [0] + list(range(1_000_000, 1_000_000 + 100 * BLOCK, 100))
    ends = [2_300_000] + [s + 500 for s in starts[1:]]
    s, e, bm = build(starts, ends)
    got = _overlaps(s, e, bm, 2_000_000, 2_000_100).tolist()
    assert got == brute(starts, ends, 2_000_000, 2_000_100)
    assert 0 in got


@pytest.mark.parametrize("seed", range(6))
def test_matches_brute_force(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(1, 4 * BLOCK))
    starts = np.sort(rng.integers(0, 1_000_000, n)).astype(np.int64)
    lengths = rng.integers(1, 5_000, n).astype(np.int64)
    if n > 40:                                  # une pincée de très longs
        long_ix = rng.choice(n, n // 20, replace=False)
        lengths[long_ix] = rng.integers(50_000, 400_000, long_ix.size)
    ends = starts + lengths
    bm = _blockmax(ends)

    sl, el = starts.tolist(), ends.tolist()
    for _ in range(60):
        qs = int(rng.integers(0, 1_000_000))
        qe = qs + int(rng.integers(1, 60_000))
        assert _overlaps(starts, ends, bm, qs, qe).tolist() == brute(sl, el, qs, qe)


def test_results_are_sorted_by_start():
    rng = np.random.default_rng(11)
    starts = np.sort(rng.integers(0, 100_000, 2_000)).astype(np.int64)
    ends = starts + rng.integers(1, 20_000, 2_000)
    bm = _blockmax(ends)
    idx = _overlaps(starts, ends, bm, 40_000, 60_000)
    assert np.all(np.diff(idx) > 0)


def test_roundtrip_on_disk(tmp_path):
    rows = [
        ("chr1", 100, 200, {"name": "a"}),
        ("chr1", 150, 400, {"name": "b"}),
        ("chr2", 100, 200, {"name": "c"}),
        ("chr1", 900, 950, {"name": "d"}),
    ]
    TrackWriter("t", tmp_path).write(rows)
    track = Track(tmp_path / "t")
    try:
        assert len(track) == 4
        assert track.chroms == ["chr1", "chr2"]
        assert [f.attrs["name"] for f in track.query("chr1", 0, 1000)] == ["a", "b", "d"]
        assert [f.attrs["name"] for f in track.query("chr2", 0, 1000)] == ["c"]
        assert track.query("chr3", 0, 1000) == []       # chromosome absent
        assert track.count("chr1", 160, 170) == 2
    finally:
        track.close()


def test_writer_rejects_empty_interval(tmp_path):
    with pytest.raises(ValueError, match="vide ou inversé"):
        TrackWriter("t", tmp_path).write([("chr1", 500, 500, {})])


def test_chroms_are_naturally_ordered(tmp_path):
    rows = [(c, 10, 20, {}) for c in ("chr10", "chrX", "chr2", "chrM", "chr1", "chrY")]
    TrackWriter("t", tmp_path).write(rows)
    track = Track(tmp_path / "t")
    try:
        assert track.chroms == ["chr1", "chr2", "chr10", "chrX", "chrY", "chrM"]
    finally:
        track.close()


def test_missing_store_says_what_to_do(tmp_path):
    with pytest.raises(FileNotFoundError, match="geno build"):
        Store(tmp_path / "nowhere")


def test_coalesced_reads_match_single_reads(tmp_path):
    """La lecture groupée doit rendre exactement ce que rendent les lectures une à une.

    C'est l'invariant que l'optimisation ne doit pas casser : mêmes attributs,
    même ordre, quels que soient les trous entre enregistrements.
    """
    rows = [("chr1", i * 1000, i * 1000 + 500, {"i": i, "pad": "x" * (i % 97)}) for i in range(5000)]
    TrackWriter("t", tmp_path).write(rows)
    track = Track(tmp_path / "t")
    try:
        for lo, hi in ((0, 5_000_000), (1_200_000, 1_260_000), (4_999_000, 5_000_000)):
            feats = track.query("chr1", lo, hi)
            expected = [track._attrs(f.attrs["i"]) for f in feats]
            assert [f.attrs for f in feats] == expected
            assert [f.attrs["i"] for f in feats] == sorted(f.attrs["i"] for f in feats)
    finally:
        track.close()


def test_sparse_hits_do_not_read_the_whole_blob(tmp_path):
    """Deux hits éloignés dans l'index ne doivent pas déclencher la lecture de tout le blob.

    Le cas réel : un très long intervalle en tête de chromosome plus un hit local
    très loin derrière. Leurs enregistrements sont séparés par tous ceux du
    milieu, qu'il ne faut surtout pas lire.
    """
    rows = [("chr1", 0, 900_000_000, {"pad": "a" * 200})]                       # index 0
    rows += [("chr1", i * 1000 + 1000, i * 1000 + 1050, {"pad": "b" * 200})     # le milieu
             for i in range(5000)]
    rows.append(("chr1", 899_999_000, 900_000_000, {"pad": "c" * 200}))         # dernier index
    TrackWriter("t", tmp_path).write(rows)
    track = Track(tmp_path / "t")
    try:
        reads = []
        real = track._blob.read
        track._blob.read = lambda n: (reads.append(n), real(n))[1]
        feats = track.query("chr1", 899_999_500, 899_999_600)
        blob_size = (tmp_path / "t" / "attrs.jsonl").stat().st_size
        assert len(feats) == 2                       # le long et le local
        assert len(reads) == 2                       # une lecture chacun, pas une seule énorme
        assert sum(reads) < blob_size // 100         # on n'a pas traversé le milieu
    finally:
        track._blob.read = real
        track.close()


def test_contiguous_hits_are_read_in_one_go(tmp_path):
    """Le pendant : des hits contigus ne doivent PAS coûter un aller-retour chacun."""
    rows = [("chr1", i * 100, i * 100 + 150, {"pad": "x" * 50}) for i in range(2000)]
    TrackWriter("t", tmp_path).write(rows)
    track = Track(tmp_path / "t")
    try:
        reads = []
        real = track._blob.read
        track._blob.read = lambda n: (reads.append(n), real(n))[1]
        feats = track.query("chr1", 0, 200_000)
        assert len(feats) == 2000
        assert len(reads) == 1
    finally:
        track._blob.read = real
        track.close()
