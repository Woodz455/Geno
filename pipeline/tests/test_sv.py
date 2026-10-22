"""Mode variant structural (semaine 15) : ce qui se teste sans les sources réelles.

Le locus réel et ses allèles 3D sont vérifiés par `geno sv check` ; ici, les briques contre
des oracles écrits à la main.
"""

from __future__ import annotations

from collections import deque

import numpy as np

from geno_pipeline.fine.region import Region
from geno_pipeline.sv import locus as L
from geno_pipeline.sv import predict as P


def _bfs(n, loops, s):
    adj = [[] for _ in range(n)]
    for i in range(n - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    for i, j in loops:
        adj[i].append(j)
        adj[j].append(i)
    d = [-1] * n
    d[s] = 0
    q = deque([s])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if d[v] < 0:
                d[v] = d[u] + 1
                q.append(v)
    return np.array(d, float)


def test_a_held_loop_is_a_shortcut():
    legs = np.array([[10, 50], [-1, -1]])
    d = P.distances(100, legs[None], np.array([0]))[0, 0]
    assert d[60] == 10 + 1 + 10 and d[50] == 11 and d[30] == 30


def test_distances_match_brute_force_on_nested_and_disjoint_loops():
    rng = np.random.default_rng(1)
    n = 300
    loops = [(20, 120), (40, 90), (50, 60), (150, 260), (200, 210)]
    legs = np.array(loops + [(-1, -1)] * 3)
    src = rng.choice(n, 12, replace=False)
    d = P.distances(n, legs[None], src)[0]
    for k, s in enumerate(src):
        assert np.array_equal(d[k], _bfs(n, loops, int(s)))


def test_law_is_capped_and_decreasing():
    d = np.array([0.0, 0.5, 1.0, 2.0, 8.0])
    p = P.law(d, 1.5, 1.0)
    assert p[0] == 1 and p[1] == 1 and p[2] == 1 and p[3] < 1 and p[4] < p[3]


def test_calibration_recovers_a_planted_law():
    # Des « observations » tirées de la loi elle-même : le calage doit la retrouver.
    rng = np.random.default_rng(2)
    n, frames = 200, 40
    legs = np.full((frames, 4, 2), -1)
    for f in range(frames):
        i = rng.integers(0, n - 60)
        legs[f, 0] = (i, i + rng.integers(10, 50))
    src = np.arange(0, n, 10)
    truth = P.rows(n, legs, src, 1.6, 1.4)
    coords = np.zeros((1, n, 3))
    obs_orig = P.observed
    try:
        P.observed = lambda c, s, cutoff=1.5: truth
        fit = P.calibrate(coords.repeat(frames, axis=0), legs, every=10, frames=frames)
    finally:
        P.observed = obs_orig
    assert abs(fit.alpha - 1.6) <= 0.05 and abs(fit.d0 - 1.4) / 1.4 < 0.1 and fit.r_log > 0.999


def _locus():
    region = Region("chr2", 0, 4_000_000, 2_000)
    genes = {
        "EPHA4": L.Gene("EPHA4", 2_500_000, 2_600_000, "-"),     # TSS 2 599 999
        "PAX3": L.Gene("PAX3", 3_300_000, 3_400_000, "-"),
    }
    peaks = [L.Peak(700_000, 700_400, "+", 300.0), L.Peak(2_600_100, 2_600_500, "-", 274.0),
             L.Peak(2_601_000, 2_601_300, "-", 30.0), L.Peak(3_500_000, 3_500_400, "+", 100.0),
             L.Peak(1_500_000, 1_500_400, None, 500.0)]
    loc = L.Locus(region, genes, peaks, 274.0, (920_000, 950_000), (2_600_100, 2_601_300), (700_000, 2_601_300))
    loc.variants = L.family(loc)
    return loc


def test_every_del_b_obeys_the_text():
    loc = _locus()
    dels = [v for v in loc.variants if v.kind == "DEL-B"]
    assert dels
    for v in dels:
        assert v.length in L.DEL_LENGTHS
        assert v.start >= loc.enhancer[1]                        # les enhancers survivent
        assert loc.boundary[1] < v.end < loc.genes["PAX3"].start   # la frontière part, PAX3 reste
    for v in loc.variants:
        if v.kind == "DelBS":
            assert v.end < loc.boundary[0]                          # la frontière reste


def test_deletion_remaps_sites_and_the_spared_boundary_stays_at_the_junction():
    loc = _locus()
    v = next(x for x in loc.variants if x.kind == "DEL-B")
    keep = next(x for x in loc.variants if x.kind == "DEL-B, frontière épargnée" and x.start == v.start and x.end == v.end)
    s0, m0 = loc.sites(None)
    s1, m1 = loc.sites(v)
    s2, m2 = loc.sites(keep)
    k = int(((m0 >= v.start) & (m0 < v.end)).sum())
    assert len(m1) == len(m0) - k == len(m2)
    assert s0.k == 4                                          # le pic sans motif n'est pas un site
    assert s1.k == 2                                          # la frontière est partie avec la délétion
    assert s2.k == 4 and set(s2.pos) >= set(s1.pos)
    junction = int(np.flatnonzero(m2 < v.start).max())
    assert sorted(s2.pos)[1] == junction == sorted(s2.pos)[2]
    # Un site après la délétion recule d'exactement k monomères.
    assert max(s1.pos) == max(s0.pos) - k
    assert loc.allele_region(v).n == len(m1)


def test_occupancy_follows_the_signal_with_floor_and_ceiling():
    loc = _locus()
    assert loc.occupancy(274.0) == 0.95 and loc.occupancy(1_000.0) == 0.95
    assert abs(loc.occupancy(137.0) - 0.475) < 1e-9 and loc.occupancy(1.0) == 0.05
