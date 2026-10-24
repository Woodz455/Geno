"""Les vues liées (semaine 14) : la carte de contacts sur la grille de HiGlass.

Le fichier réel est vérifié par `geno linked check` (conservation, pixels recomposés par
distances directes, témoin décalé) ; ici, les briques, chacune contre une force brute.
"""

from __future__ import annotations

import random

import numpy as np
from scipy import sparse

from geno_pipeline.export import linked as L


def _brute_spread(starts, ends, res, n_px):
    out = np.zeros((len(starts), n_px))
    for k, (s, e) in enumerate(zip(starts, ends)):
        for p in range(n_px):
            out[k, p] = max(0, min(e, (p + 1) * res) - max(s, p * res)) / (e - s)
    return out


def test_spread_matches_brute_force_and_rows_sum_to_one():
    rng = random.Random(1)
    starts, ends, at = [], [], 0
    for _ in range(40):
        n = rng.randrange(300, 2_500)
        starts.append(at)
        ends.append(at + n)
        at += n
    n_px = -(-at // 1_000)
    s = L.spread(np.array(starts), np.array(ends), 1_000, n_px).toarray()
    assert np.allclose(s, _brute_spread(starts, ends, 1_000, n_px))
    assert np.allclose(s.sum(axis=1), 1.0)


def test_pixels_conserve_contacts_and_match_brute_force():
    rng = np.random.default_rng(2)
    n = 30
    starts = np.cumsum(rng.integers(500, 1_500, n)) - 500
    ends = starts + rng.integers(400, 900, n)
    m = rng.random((n, n))
    m = (m + m.T) / 2
    n_px = -(-int(ends.max()) // 700)
    s = L.spread(starts, ends, 700, n_px)
    px = L.to_pixels(sparse.csr_matrix(m), s).toarray()
    b = _brute_spread(starts, ends, 700, n_px)
    assert np.allclose(px, b.T @ m @ b)
    assert np.isclose(px.sum(), m.sum())


def test_symmetric_counts_the_diagonal_once():
    m = L.symmetric(np.array([0, 0, 1]), np.array([0, 2, 1]), np.array([5.0, 2.0, 1.0]), 3).toarray()
    assert m.tolist() == [[5, 0, 2], [0, 1, 0], [2, 0, 0]]


def test_tile_layout_is_column_x_row_y_and_only_lower_half_is_written():
    n_px = 600                              # 3 tuiles de côté
    rng = np.random.default_rng(3)
    a = rng.random((n_px, n_px))
    a = sparse.csr_matrix(np.tril(a) + np.tril(a, -1).T)
    got = {(x, y): d for x, y, d in L.tiles(a, n_px)}
    assert all(x <= y for x, y in got)
    assert set(got) == {(x, y) for y in range(3) for x in range(y + 1)}
    full = a.toarray()
    for (x, y), d in got.items():
        blk = full[y * 256:(y + 1) * 256, x * 256:(x + 1) * 256]
        assert np.allclose(d[:blk.shape[0], :blk.shape[1]], blk)
        assert (d[blk.shape[0]:] == 0).all() and (d[:, blk.shape[1]:] == 0).all()


def test_region_masks_leave_cross_region_pixels_nan():
    n_px = 300
    m1 = np.zeros(n_px, bool)
    m1[10:20] = True
    m2 = np.zeros(n_px, bool)
    m2[200:210] = True
    a = sparse.csr_matrix(np.ones((n_px, n_px)))
    got = {(x, y): d for x, y, d in L.tiles(a, n_px, [m1, m2])}
    d00 = got[(0, 0)]
    assert d00[15, 12] == 1 and np.isnan(d00[15, 100]) and np.isnan(d00[205, 15])
    assert d00[205, 202] == 1


def test_tile_round_trip_is_exact():
    rng = np.random.default_rng(4)
    d = (rng.random((256, 256)) * 300).astype(np.float32)
    d[3, 4] = np.nan
    back = L.decode_tile(L.encode_tile(d))
    assert np.array_equal(np.isnan(back), np.isnan(d))
    ok = ~np.isnan(d)
    assert np.allclose(back[ok], d[ok], rtol=4.9e-4, atol=0)        # float16 : 11 bits


def test_tileset_totals_are_conserved_at_every_resolution():
    starts = np.arange(0, 50_000, 1_000)
    ends = starts + 1_000
    rng = np.random.default_rng(5)
    v = rng.random((50, 50))
    src = L.Source(starts, ends, sparse.csr_matrix((v + v.T) / 2), 10, 1.5, "test")
    blobs = []

    def put(b):
        blobs.append(b)
        return {"i": len(blobs) - 1}

    ts = L.tileset("t", "test", [src], (8_000, 3_000, 1_000), 50_000, put)
    for lv in ts["levels"]:
        s = 0.0
        for key, spec in lv["tiles"].items():
            _, x, y = map(int, key.split("."))
            t = np.nan_to_num(L.decode_tile(blobs[spec["i"]]).astype(np.float64))
            s += t.sum() * (1 if x == y else 2)
        assert np.isclose(s, src.matrix.sum(), rtol=1e-3)
