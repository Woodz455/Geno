"""Semaine 21 — les outils de validation, sur des données plantées (aucun fichier réel requis)."""

import numpy as np

from geno_pipeline.validate.common import (Model, null_p95, over_expected, pair_distances, spearman,
                                           verdict)
from geno_pipeline.validate.sprite import enrichment, oe_vec


def toy_model(seed=0, structures=20, beads=40):
    """Deux copies d'un chromosome, marches aléatoires dans un noyau de rayon 5."""
    rng = np.random.default_rng(seed)
    n = 2 * beads
    steps = rng.normal(scale=0.2, size=(structures, n, 3))
    coords = np.cumsum(steps, axis=1).astype(np.float32)
    start = np.tile(np.arange(beads) * 750_000, 2)
    return Model(coords=coords, radius=np.full(n, 0.17), copy_id=np.repeat([0, 1], beads),
                 start=start, end=start + 750_000, labels=("chr21:a", "chr21:b"), nuclear_radius=5.0,
                 lad_source="synthetic", lad=rng.random(n))


def test_verdicts():
    assert verdict(0.5, 0.3, p95=0.1) == "accord"
    assert verdict(0.2, 0.3, p95=0.1) == "partiel"
    assert verdict(0.2, 0.3, p95=0.25) == "échec"
    assert verdict(-0.1, 0.3, p95=-0.2) == "échec"
    assert verdict(-0.5, -0.3, sign=-1, p95=0.1) == "accord"
    assert verdict(float("nan"), 0.3) == "échec"


def test_oe_and_expected():
    k = np.array([1, 1, 2, 2])
    v = np.array([1.0, 3.0, 4.0, 4.0])
    np.testing.assert_allclose(oe_vec(v, k), [0.5, 1.5, 1.0, 1.0])
    m = np.array([[np.nan, 2.0], [2.0, np.nan]])
    off = np.array([[0, 1], [1, 0]])
    np.testing.assert_allclose(over_expected(m, off)[0, 1], 1.0)


def test_planted_data_equal_to_model_is_accord_and_permuted_is_not():
    model = toy_model()
    idx = model.beads_at("chr21", (np.arange(40) * 750_000 + 375_000))
    med = np.nanmedian(pair_distances(model, idx, cis_only=True), axis=0)
    i, j = np.triu_indices(40, 1)
    k = (j - i).astype(np.int64)
    data = oe_vec(med[i, j], k)                         # les « données » sont le modèle lui-même
    r = spearman(oe_vec(med[i, j], k), data)
    p95, mean = null_p95(lambda v: spearman(v, data),
                         lambda rng: (lambda p: oe_vec(med[np.ix_(p, p)][i, j], k))(rng.permutation(40)), n=200)
    assert r > 0.999 and verdict(r, 0.3, p95=p95) == "accord"
    assert abs(mean) < 0.1                              # les permutations sont centrées sur 0
    rng = np.random.default_rng(1)
    shuffled = data[rng.permutation(len(data))]
    assert verdict(spearman(oe_vec(med[i, j], k), shuffled), 0.3, p95=p95) != "accord"


def test_beads_at_and_cis_pairs():
    model = toy_model()
    idx = model.beads_at("chr21", np.array([10, 750_010, 29_999_999, 40_000_000]))
    assert idx.tolist() == [[0, 1, 39, -1], [40, 41, 79, -1]]
    d = pair_distances(model, idx[:, :2], cis_only=True)
    assert d.shape == (40, 2, 2)                        # 2 copies × 20 repliements
    assert np.all(np.isfinite(d))


def test_imposed_radial_correlation_is_seen():
    # Témoin positif de C5 : un rayon qui suit la piste doit donner ρ ≈ −1.
    lad = np.linspace(0, 1, 50)
    radial = 1.0 - lad + np.random.default_rng(2).normal(scale=0.01, size=50)
    assert spearman(radial, lad) < -0.95


def test_enrichment_finds_a_planted_pair():
    R = np.array([3.0, 2.0, 1.5, 1.0])
    M = np.outer(R, R)
    M[0, 2] = M[2, 0] = 3 * M[0, 2]                       # une paire enrichie exprès
    E = enrichment(M)
    iu = np.triu_indices(4, 1)
    assert (iu[0][np.argmax(E[iu])], iu[1][np.argmax(E[iu])]) == (0, 2)
