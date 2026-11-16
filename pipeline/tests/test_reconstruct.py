"""Contact → distance → 3D, validé contre une conformation connue.

La semaine 3 plantait des enrichissements de contacts ; ici on plante une vraie
géométrie, on en dérive les contacts, et on vérifie qu'on la retrouve. C'est ce
que les données réelles ne permettront jamais : dans du Hi-C réel, la structure
3D est précisément l'inconnue.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy", reason="pile Hi-C absente — voir docs/SETUP.md")

from geno_pipeline.hic.polymer import chain, contacts  # noqa: E402
from geno_pipeline.hic.reconstruct import (  # noqa: E402
    balance_dense,
    best,
    classical_mds,
    complete_shortest_path,
    fidelity,
    procrustes,
    reconstruct,
    sweep,
    to_distance,
)

N = 180
ALPHAS = np.round(np.arange(0.20, 0.70, 0.05), 3)


@pytest.fixture(scope="module")
def conf():
    return chain(n=N, seed=1)


# --------------------------------------------------------------------------
# La conformation plantée
# --------------------------------------------------------------------------


def test_chain_is_a_chain_not_a_cloud(conf):
    bonds = np.linalg.norm(np.diff(conf.coords, axis=0), axis=1)
    assert 0.5 < bonds.mean() < 1.6, f"pas de liaison moyenne {bonds.mean():.2f}"
    assert bonds.max() < 3.0, "une liaison a cassé"


def test_confinement_is_respected(conf):
    assert np.linalg.norm(conf.coords, axis=1).max() <= conf.radius * 1.02


def test_beads_do_not_interpenetrate(conf):
    d = conf.distances()
    np.fill_diagonal(d, np.inf)
    # Les voisins de chaîne se touchent par construction ; on regarde les autres.
    far = np.abs(np.arange(N)[:, None] - np.arange(N)[None, :]) > 2
    assert d[far].min() > 0.3, "des billes non voisines se superposent"


@pytest.mark.parametrize("n,seed", [(150, 0), (180, 3), (350, 1), (500, 2)])
def test_compartment_a_is_central(n, seed):
    """Territoires de Cremer & Cremer : A au centre, B en périphérie.

    C'est le fait indépendant que la structure plantée doit reproduire, et que la
    reconstruction devra retrouver.

    Testé sur plusieurs tailles et graines délibérément. La première version de
    `chain()` obtenait la ségrégation par une dérive pendant la marche, qui luttait
    contre la diffusion du hasard : l'écart valait +0,157 à 180 billes et +0,009 à
    350. Un test à une seule taille l'aurait déclarée acquise.
    """
    conf = chain(n=n, seed=seed)
    r = conf.radial()
    gap = r[conf.compartment == -1].mean() - r[conf.compartment == 1].mean()
    assert gap > 0.10, f"ségrégation A/B trop faible à n={n}, graine={seed} : {gap:+.3f}"


# --------------------------------------------------------------------------
# Les briques, séparément
# --------------------------------------------------------------------------


def test_zero_contacts_become_infinite_not_zero():
    counts = np.array([[0.0, 4.0, 0.0], [4.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    d = to_distance(counts, 0.5)
    assert np.isinf(d[0, 2]), "une paire sans contact n'est pas à distance nulle"
    assert d[0, 1] < d[1, 2], "plus de contacts doit vouloir dire plus proche"
    assert (np.diag(d) == 0).all()


def test_shortest_path_completion_removes_infinities_and_respects_triangle():
    counts = np.array([[0.0, 9.0, 0.0], [9.0, 0.0, 9.0], [0.0, 9.0, 0.0]])
    geo = complete_shortest_path(to_distance(counts, 0.5))
    assert np.isfinite(geo).all()
    assert geo[0, 2] == pytest.approx(geo[0, 1] + geo[1, 2])
    n = len(geo)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                assert geo[i, j] <= geo[i, k] + geo[k, j] + 1e-9


def test_mds_recovers_exact_euclidean_coordinates():
    rng = np.random.default_rng(0)
    truth = rng.normal(0, 1, (40, 3))
    d = np.linalg.norm(truth[:, None] - truth[None, :], axis=2)
    _aligned, rmsd = procrustes(classical_mds(d, 3), truth)
    assert rmsd < 1e-8, f"MDS exact devrait être exact (rmsd={rmsd:.2e})"


def test_procrustes_is_blind_to_rotation_translation_and_scale():
    rng = np.random.default_rng(3)
    truth = rng.normal(0, 1, (30, 3))
    q, _ = np.linalg.qr(rng.normal(0, 1, (3, 3)))
    moved = 7.5 * truth @ q + np.array([100.0, -3.0, 42.0])
    _a, rmsd = procrustes(moved, truth)
    assert rmsd < 1e-8


def test_procrustes_accepts_a_mirror_image():
    """Une matrice de distances ne porte aucune chiralité.

    Une reconstruction miroir est donc une reconstruction correcte. Interdire la
    réflexion compterait la moitié des solutions valides comme des échecs.
    """
    rng = np.random.default_rng(4)
    truth = rng.normal(0, 1, (30, 3))
    mirrored = truth * np.array([1.0, 1.0, -1.0])
    _a, rmsd = procrustes(mirrored, truth)
    assert rmsd < 1e-8


def test_balancing_flattens_dense_marginals():
    rng = np.random.default_rng(5)
    base = rng.random((60, 60)) + 0.5
    base = (base + base.T) / 2
    np.fill_diagonal(base, 0)
    bias = np.exp(rng.normal(0, 0.7, 60))
    skewed = base * np.outer(bias, bias)
    cv = lambda m: float(np.std(m.sum(1)) / np.mean(m.sum(1)))  # noqa: E731
    assert cv(balance_dense(skewed)) < cv(skewed) / 5


# --------------------------------------------------------------------------
# La chaîne complète
# --------------------------------------------------------------------------


def test_reconstruction_recovers_the_planted_geometry(conf):
    counts = contacts(conf, gamma=3.0, total=4_000_000, seed=2)
    f = best(sweep(counts, conf.coords, ALPHAS))
    assert f.nrmsd < 0.30, str(f)
    assert f.dist_rho > 0.95, str(f)
    assert f.radial_r > 0.90, "l'organisation radiale doit survivre"


def test_reconstruction_keeps_a_central_and_b_peripheral(conf):
    """Le fait indépendant, retrouvé de l'autre côté de la reconstruction."""
    counts = contacts(conf, gamma=3.0, total=8_000_000, seed=7)
    coords = reconstruct(counts, 0.40)
    centre = coords.mean(axis=0)
    r = np.linalg.norm(coords - centre, axis=1)
    assert r[conf.compartment == 1].mean() < r[conf.compartment == -1].mean()


def test_geodesic_completion_is_not_optional(conf):
    """La complétion par plus courts chemins n'est pas un raffinement.

    Remplacer les paires sans contact par une grande constante au lieu de leur
    géodésique dégrade la reconstruction d'un facteur net. C'est l'apport de
    ShRec3D, et il se mesure.
    """
    counts = contacts(conf, gamma=3.0, total=1_500_000, seed=2)
    alpha = 0.45

    with_geo = fidelity(reconstruct(counts, alpha), conf.coords, alpha)

    d = to_distance(counts, alpha)
    capped = np.where(np.isfinite(d), d, d[np.isfinite(d)].max() * 2.0)
    np.fill_diagonal(capped, 0.0)
    without = fidelity(classical_mds(capped), conf.coords, alpha)

    assert with_geo.nrmsd < without.nrmsd / 1.4, f"avec {with_geo.nrmsd:.3f}, sans {without.nrmsd:.3f}"


def test_reconstruction_improves_with_depth(conf):
    """Plus de contacts, meilleure géométrie. Monotone, et c'est la seule chose
    qui le soit vraiment dans ce balayage."""
    errs = []
    for total in (500_000, 5_000_000, 50_000_000):
        counts = contacts(conf, gamma=3.0, total=total, seed=2)
        errs.append(best(sweep(counts, conf.coords, ALPHAS)).nrmsd)
    assert errs[0] > errs[1] > errs[2], f"nRMSD doit décroître avec la profondeur : {errs}"


def test_alpha_differs_from_one_over_gamma_at_finite_depth(conf):
    """À profondeur finie, l'exposant optimal n'est pas 1/gamma.

    **Attention à ce que ce test n'affirme pas.** Une première version affirmait
    que l'écart est positif et décroît avec la profondeur. C'était trop fort : en
    changeant la longueur de chaîne, l'optimum à faible profondeur passait de
    0,250 à 0,575. Le minimum est plat quand les données sont creuses, et l'argmin
    y devient instable.

    Ce qui tient : l'optimum s'écarte de 1/gamma, et reprendre alpha = 1/3 d'un
    article sans le calibrer sur ses propres données est injustifié.
    """
    counts = contacts(conf, gamma=3.0, total=1_500_000, seed=2)
    f = best(sweep(counts, conf.coords, ALPHAS))
    assert abs(f.alpha - 1 / 3) > 0.04, f"attendu un écart mesurable à 0,333 : {f}"


def test_noiseless_limit_lands_on_one_over_gamma(conf):
    """Sans bruit, l'inversion exacte doit gagner — et c'est bien 1/gamma."""
    d = conf.distances()
    np.fill_diagonal(d, np.inf)
    exact = d ** (-3.0)
    np.fill_diagonal(exact, 0.0)
    fine = np.round(np.arange(0.25, 0.50, 0.025), 3)
    f = best([fidelity(reconstruct(exact, a), conf.coords, a) for a in fine])
    assert abs(f.alpha - 1 / 3) <= 0.03, str(f)
    assert f.nrmsd < 0.02, f"sans bruit, la reconstruction doit être quasi exacte : {f}"
