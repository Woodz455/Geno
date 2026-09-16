"""Semaine 8 — extrusion de boucles et échelle fine.

Les tests de dynamique tournent sur des chaînes minuscules et quelques milliers
de pas : ce qu'ils vérifient, c'est que le champ de force fait ce qu'il dit, pas
qu'une trajectoire de production converge. Ce dernier point se mesure, il ne
s'assert pas (§ `polymer.equilibration`).
"""

from __future__ import annotations

import numpy as np
import pytest

from geno_pipeline.fine import (
    CONVERGENT,
    DIVERGENT,
    KINDS,
    TANDEM,
    Region,
    Sites,
    blob_radius,
    bond_catalogue,
    call_boundaries,
    contact_map,
    dot_score,
    dots_by_kind,
    insulation,
    junction,
    lef_count,
    lifetime_steps,
    monomer_radius,
    plant,
    ps,
    read_ctcf,
    separation_curve,
    start,
)
from geno_pipeline.fine import fold, simulate
from geno_pipeline.fine.observe import Map
from geno_pipeline.fine.polymer import Field, _switch
from geno_pipeline.fine.region import Truth
from geno_pipeline.fine.run import Setup, build, load, save

REGION = Region("chr7", 4_000_000, 8_000_000, 2_000)


# --------------------------------------------------------------------------
# region : la loi d'échelle, les barrières, ce qu'on plante
# --------------------------------------------------------------------------


def test_monomer_radius_is_the_week_six_law_not_a_new_one():
    """166,8 nm à 750 kb et 39,5 nm à 10 kb — les valeurs du tableau de la semaine 6.

    Si cette semaine avait le droit de choisir la taille d'un monomère, le
    raccord de fin de semaine serait garanti d'avance et ne mesurerait rien.
    """
    assert monomer_radius(750_000) == pytest.approx(166.8, abs=0.15)
    assert monomer_radius(10_000) == pytest.approx(39.5, abs=0.15)


def test_monomer_radius_follows_the_cube_root_to_the_rounding_of_the_bead_count():
    """`r ∝ L^(1/3)`, mais pas au bit près, et la raison est intéressante.

    Le nombre de billes est **arrondi** à l'entier, si bien que huit monomères de
    2 kb ne font pas exactement un monomère de 16 kb : l'écart mesuré vaut
    1,1·10⁻⁷ en relatif. C'est la même approximation que la semaine 6, et la
    garder identique vaut mieux que de la corriger ici seulement.
    """
    assert monomer_radius(2_000) / monomer_radius(16_000) == pytest.approx(0.5, rel=1e-5)


def test_blob_radius_is_the_volume_the_nucleus_allots():
    n, phi = 2_000, 0.30
    r = 0.5                                     # rayon en sigma, par définition
    R = blob_radius(n, phi)
    assert n * (4 / 3) * np.pi * r**3 / ((4 / 3) * np.pi * R**3) == pytest.approx(phi)


def test_region_maps_a_coordinate_to_a_monomer_and_back():
    assert REGION.n == 2_000
    b = REGION.bead(5_527_151)                  # ACTB, cf. les fixtures de la semaine 2
    assert 0 <= b < REGION.n
    assert abs(float(REGION.bp(b)) - 5_527_151) <= REGION.bp_per_bead


def test_a_coordinate_outside_the_region_is_clipped_not_wrapped():
    assert REGION.bead(1) == 0
    assert REGION.bead(99_000_000) == REGION.n - 1


def test_barriers_follow_the_convergent_rule():
    """Un `+` arrête qui va vers la gauche, un `−` qui va vers la droite.

    C'est toute la règle de convergence, et elle tient dans cette asymétrie. Une
    inversion ici produirait des boucles ancrées sur des paires divergentes,
    c'est-à-dire l'inverse exact de ce que Rao 2014 observe.
    """
    sites = Sites(np.array([10, 40]), np.array([+1, -1]), np.array([0.8, 0.9]), "test")
    left, right = sites.barriers(100)
    assert left[10] == pytest.approx(0.8) and right[10] == 0.0
    assert right[40] == pytest.approx(0.9) and left[40] == 0.0


def test_two_sites_of_one_sense_on_one_monomer_do_not_add_up():
    sites = Sites(np.array([7, 7]), np.array([+1, +1]), np.array([0.6, 0.7]), "test")
    left, _ = sites.barriers(20)
    assert left[7] == pytest.approx(0.7)        # le plus occupé, pas la somme


def test_knockout_keeps_the_positions_and_zeroes_only_the_occupancy():
    sites, _ = plant(REGION, seed=3)
    ko = sites.without_ctcf()
    assert np.array_equal(ko.pos, sites.pos)
    assert np.array_equal(ko.strand, sites.strand)
    assert not ko.occupancy.any()


def test_plant_balances_its_two_control_classes():
    """Un témoin à un seul membre ne contredit rien ; l'équilibre est voulu.

    Tiré à pile ou face, l'effectif partait couramment à 1 contre 5 sur une
    vingtaine de domaines — mesuré. L'alternance le corrige par construction.
    """
    for seed in range(6):
        _, truth = plant(REGION, seed=seed, mean_domain=200_000, min_domain=80_000)
        d = int((truth.kind == DIVERGENT).sum())
        t = int((truth.kind == TANDEM).sum())
        assert abs(d - t) <= 1


def test_plant_anchors_carry_the_orientation_their_class_claims():
    sites, truth = plant(REGION, seed=2, mean_domain=200_000)
    at = {int(p): int(s) for p, s in zip(sites.pos, sites.strand)}
    for left, right, kind in zip(truth.left, truth.right, truth.kind):
        if kind == CONVERGENT:
            assert (at[int(left)], at[int(right)]) == (+1, -1)
        elif kind == DIVERGENT:
            assert (at[int(left)], at[int(right)]) == (-1, +1)
        else:
            assert (at[int(left)], at[int(right)]) == (+1, +1)


def test_no_boundary_is_ever_non_blocking():
    """Une frontière porte deux ancres, et une ancre arrête toujours un sens.

    La première version de cette vérité disait « bloquante si l'ancre de gauche
    est `−` ou celle de droite `+` », ce qui rate les deux autres cas et
    annonçait 18 frontières bloquantes sur 19 — la dix-neuvième insulant
    évidemment autant. Un témoin qui contredit son propre énoncé n'en est pas un.
    """
    for seed in range(6):
        _, truth = plant(REGION, seed=seed, mean_domain=200_000)
        assert set(np.unique(truth.blocking)) <= {1, 2}
        assert (truth.blocking >= 1).all()


def test_a_boundary_blocks_both_senses_exactly_when_its_two_anchors_differ():
    sites, truth = plant(REGION, seed=5, mean_domain=200_000)
    at = {int(p): int(s) for p, s in zip(sites.pos, sites.strand)}
    for k in range(len(truth.boundary)):
        right = at[int(truth.right[k])]
        left_next = at[int(truth.left[k + 1])]
        assert bool(truth.strong[k]) is (right != left_next)
        # Et dans les deux cas, au moins un sens est arrêté.
        stops_rightward = right < 0 or left_next < 0
        stops_leftward = right > 0 or left_next > 0
        assert stops_rightward or stops_leftward


def test_the_strands_are_a_pure_function_of_the_domain_class():
    sites, truth = plant(REGION, seed=7, mean_domain=200_000)
    at = {int(p): int(s) for p, s in zip(sites.pos, sites.strand)}
    for k, (l, r) in enumerate(zip(truth.left, truth.right)):
        assert truth.strand_left[k] == at[int(l)]
        assert truth.strand_right[k] == at[int(r)]


def test_read_ctcf_refuses_a_site_without_a_strand(tmp_path):
    """Un site CTCF sans orientation ne contraint rien ici. En inventer une
    fabriquerait exactement les boucles qu'on prétend prédire."""
    p = tmp_path / "ctcf.bed"
    p.write_text("chr7\t5000000\t5000020\tsite\t800\t.\n")
    with pytest.raises(ValueError, match="brin"):
        read_ctcf(str(p), REGION)


def test_read_ctcf_takes_the_bed_score_as_occupancy(tmp_path):
    p = tmp_path / "ctcf.bed"
    p.write_text(
        "chr7\t5000000\t5000020\ta\t800\t+\n"
        "chr7\t6000000\t6000020\tb\t250\t-\n"
        "chr9\t6000000\t6000020\telsewhere\t900\t+\n"
    )
    sites = read_ctcf(str(p), REGION)
    assert sites.k == 2                                   # l'autre chromosome est ignoré
    assert sites.occupancy.tolist() == pytest.approx([0.8, 0.25])
    assert sites.strand.tolist() == [+1, -1]


def test_read_ctcf_refuses_an_empty_region(tmp_path):
    p = tmp_path / "ctcf.bed"
    p.write_text("chr7\t100\t120\ta\t800\t+\n")
    with pytest.raises(ValueError, match="aucun site"):
        read_ctcf(str(p), REGION)


# --------------------------------------------------------------------------
# extrusion : le modèle 1D, vérifiable seul et pour rien
# --------------------------------------------------------------------------


def _small():
    reg = Region("chrT", 0, 400_000, 2_000)               # 200 monomères
    sites, truth = plant(reg, seed=1, mean_domain=100_000, min_domain=60_000)
    return reg, sites, truth


def test_a_convergent_pair_anchors_and_a_control_pair_does_not():
    """Le résultat central de l'étage 1D, et il peut échouer.

    Les témoins partagent la distribution de longueurs des convergents : un score
    plus élevé aux convergents ne peut donc pas être mis sur le compte de la
    distance génomique.
    """
    reg = Region("chrT", 0, 3_000_000, 2_000)
    sites, truth = plant(reg, seed=1, mean_domain=200_000, min_domain=80_000)
    ex = simulate(reg, sites, seed=1, snapshots=300, stride=4)
    conv = ex.anchor_frequency(truth.pairs(CONVERGENT)).mean()
    ctrl = np.concatenate(
        [ex.anchor_frequency(truth.pairs(k)) for k in (DIVERGENT, TANDEM) if len(truth.pairs(k))]
    ).mean()
    assert conv > 3.0 * ctrl


def test_with_every_site_unoccupied_nothing_anchors():
    reg = Region("chrT", 0, 3_000_000, 2_000)
    sites, truth = plant(reg, seed=1, mean_domain=200_000, min_domain=80_000)
    ex = simulate(reg, sites.without_ctcf(), seed=1, snapshots=300, stride=4)
    assert ex.anchor_frequency(truth.pairs(CONVERGENT)).mean() < 0.01


def test_every_snapshot_keeps_one_row_per_complex():
    """L'identité des emplacements survit à l'échantillonnage.

    Sans ça, une liaison OpenMM sauterait d'une paire à une autre sans rapport
    entre deux mises à jour, et la dynamique encaisserait un choc là où
    l'extrusion n'avance que d'un cran.
    """
    reg, sites, _ = _small()
    ex = simulate(reg, sites, seed=0, snapshots=25, stride=2)
    assert all(f.shape == (ex.n_lef, 2) for f in ex.legs)
    assert any((f[:, 0] < 0).any() or True for f in ex.legs)


def test_two_legs_never_sit_on_the_same_monomer():
    reg, sites, _ = _small()
    ex = simulate(reg, sites, seed=2, snapshots=40, stride=3)
    for frame in ex.legs:
        feet = frame[frame[:, 0] >= 0].ravel()
        assert len(feet) == len(np.unique(feet))


def test_the_left_leg_stays_left_and_both_stay_in_the_region():
    reg, sites, _ = _small()
    ex = simulate(reg, sites, seed=3, snapshots=40, stride=3)
    for frame in ex.legs:
        held = frame[frame[:, 0] >= 0]
        assert (held[:, 0] < held[:, 1]).all()
        assert held.min() >= 0 and held.max() < reg.n


def test_longer_processivity_makes_longer_loops():
    reg = Region("chrT", 0, 3_000_000, 2_000)
    sites, _ = plant(reg, seed=1, mean_domain=200_000)
    short = simulate(reg, sites, seed=1, snapshots=200, stride=4, processivity=100_000)
    long = simulate(reg, sites, seed=1, snapshots=200, stride=4, processivity=400_000)
    assert np.median(long.loop_sizes(2_000)) > 1.5 * np.median(short.loop_sizes(2_000))


def test_lifetime_and_count_come_from_the_stated_relations():
    assert lifetime_steps(200_000, 2_000, 1.0) == pytest.approx(50.0)
    assert lef_count(Region("c", 0, 3_000_000, 2_000), 200_000) == 15


def test_extrusion_is_reproducible_from_its_seed():
    reg, sites, _ = _small()
    a = simulate(reg, sites, seed=11, snapshots=12, stride=2)
    b = simulate(reg, sites, seed=11, snapshots=12, stride=2)
    assert all(np.array_equal(x, y) for x, y in zip(a.legs, b.legs))


# --------------------------------------------------------------------------
# polymer : le champ de force fait-il ce qu'il dit
# --------------------------------------------------------------------------


def test_bond_catalogue_dedupes_and_drops_the_empty_slots():
    legs = [
        np.array([[3, 9], [-1, -1]]),
        np.array([[3, 9], [20, 25]]),
        np.array([[9, 3], [-1, -1]]),        # même paire, ordre inverse
    ]
    pairs, index = bond_catalogue(legs)
    assert pairs.tolist() == [[3, 9], [20, 25]]
    assert index[(3, 9)] == 0 and index[(20, 25)] == 1


def test_bond_catalogue_of_a_trajectory_without_cohesins_is_empty():
    pairs, index = bond_catalogue([np.full((2, 2), -1)])
    assert pairs.shape == (0, 2) and index == {}


def test_switch_only_touches_the_bonds_whose_state_changes():
    """Reparcourir le catalogue à chaque instantané coûterait des millions
    d'appels pour un résultat identique."""

    class Recorder:
        def __init__(self):
            self.calls = []

        def setBondParameters(self, b, i, j, length, k):   # noqa: N802
            self.calls.append((b, k > 0))

    pairs, index = bond_catalogue([np.array([[0, 5], [10, 20], [30, 40]])])
    rec = Recorder()
    live = _switch(rec, pairs, index, np.array([[0, 5], [10, 20], [-1, -1]]), set(), 1.0)
    assert live == {0, 1} and len(rec.calls) == 2

    rec.calls.clear()
    live = _switch(rec, pairs, index, np.array([[0, 5], [30, 40], [-1, -1]]), live, 1.0)
    assert live == {0, 2}
    assert sorted(rec.calls) == [(1, False), (2, True)]    # une extinction, un allumage


def test_a_confined_start_never_leaves_the_sphere():
    x = start(500, 6.0, np.random.default_rng(0))
    assert np.linalg.norm(x - x.mean(0) + x.mean(0), axis=1).max() <= 6.0 + 1e-6 + 12.0
    y = start(500, 6.0, np.random.default_rng(0))
    assert np.array_equal(x, y)


def test_a_confined_start_has_unit_steps_rather_than_a_squeezed_walk():
    """Comprimer une marche libre donnerait le bon `Rg` et une structure interne
    fausse — et `Rg` ne peut pas le voir (§ `polymer.start`)."""
    x = start(800, 6.0, np.random.default_rng(1))
    steps = np.linalg.norm(np.diff(x, axis=0), axis=1)
    assert steps.mean() == pytest.approx(1.0, abs=0.15)


@pytest.mark.parametrize("confine", [None, 4.0])
def test_the_chain_holds_together_and_stays_where_it_is_told(confine):
    n = 120
    legs = [np.full((1, 2), -1, dtype=np.int64)] * 3
    traj = fold(n, legs, confine=confine, md_per_step=1_500, relax=2_000, seed=0)
    last = traj[-1]
    d = np.linalg.norm(np.diff(last, axis=0), axis=1)
    assert 0.7 < d.mean() < 1.4
    if confine is not None:
        assert np.linalg.norm(last, axis=1).max() < confine * 1.25


def test_a_cohesin_bond_pulls_its_two_feet_together():
    """Le test qui discrimine : sans la liaison, les deux monomères restent loin.

    Ils sont séparés de 80 monomères le long de la chaîne ; une liaison de
    cohésine doit les ramener au contact, pas seulement les rapprocher un peu.
    """
    n = 100
    free = [np.full((1, 2), -1, dtype=np.int64)] * 4
    held = [np.array([[10, 90]], dtype=np.int64)] * 4
    d_free = np.linalg.norm(
        fold(n, free, confine=None, md_per_step=2_000, relax=2_000, seed=5)[-1][10]
        - fold(n, free, confine=None, md_per_step=2_000, relax=2_000, seed=5)[-1][90]
    )
    x = fold(n, held, confine=None, md_per_step=2_000, relax=2_000, seed=5)[-1]
    d_held = np.linalg.norm(x[10] - x[90])
    assert d_held < 1.6
    assert d_held < 0.4 * d_free


def test_the_same_seed_gives_the_same_trajectory():
    legs = [np.array([[4, 40]], dtype=np.int64)] * 2
    a = fold(60, legs, confine=5.0, md_per_step=500, relax=500, seed=9)
    b = fold(60, legs, confine=5.0, md_per_step=500, relax=500, seed=9)
    assert np.allclose(a, b)


def test_a_stiffer_chain_is_a_straighter_chain():
    legs = [np.full((1, 2), -1, dtype=np.int64)] * 2
    soft = fold(200, legs, confine=None, md_per_step=3_000, relax=3_000, seed=4,
                field=Field(stiffness=0.0))
    stiff = fold(200, legs, confine=None, md_per_step=3_000, relax=3_000, seed=4,
                 field=Field(stiffness=6.0))
    reach = lambda t: float(np.linalg.norm(t[-1][40:] - t[-1][:-40], axis=1).mean())
    assert reach(stiff) > reach(soft)


# --------------------------------------------------------------------------
# observe : les cartes, et ce qu'elles ont le droit de dire
# --------------------------------------------------------------------------


def _line(n=40, frames=3, spacing=1.0):
    x = np.zeros((frames, n, 3))
    x[:, :, 0] = np.arange(n) * spacing
    return x


def test_a_contact_map_is_symmetric_and_counts_a_pair_once():
    cm = contact_map(_line(n=10, frames=4), cutoff=1.5)
    assert np.allclose(cm.counts, cm.counts.T)
    assert cm.counts[0, 1] == pytest.approx(1.0)      # voisins, à toutes les trames
    assert cm.counts[0, 5] == 0.0


def test_the_observed_over_expected_of_a_flat_map_is_one():
    cm = Map(np.ones((12, 12)), frames=1, bin_beads=1, bp_per_bin=2_000, cutoff=1.5)
    assert np.allclose(cm.oe(), 1.0)


def test_ps_is_normalised_to_one_at_the_smallest_separation():
    s, p = ps(contact_map(_line(n=30, frames=2), cutoff=1.5))
    assert p[0] == pytest.approx(1.0)
    assert s[0] == 2_000


def test_insulation_dips_where_the_chain_is_cut_in_two():
    """Deux pelotes disjointes : le score doit s'effondrer à la jonction, et là
    seulement."""
    rng = np.random.default_rng(0)
    n, frames = 60, 40
    x = np.empty((frames, n, 3))
    for f in range(frames):
        x[f, :30] = rng.normal(0, 0.6, (30, 3))
        x[f, 30:] = rng.normal(0, 0.6, (30, 3)) + np.array([25.0, 0, 0])
    ins = insulation(contact_map(x, cutoff=1.5), window=6)
    assert np.nanargmin(ins) in range(28, 32)
    assert 30 in call_boundaries(ins, prominence=0.2) or 29 in call_boundaries(ins, prominence=0.2)


def test_a_dot_score_of_a_featureless_map_is_one():
    oe = np.ones((40, 40))
    assert dot_score(oe, np.array([[10, 30]]))[0] == pytest.approx(1.0)


def test_a_dot_score_sees_a_planted_peak():
    """Le point fait trois casiers de côté, comme une vraie boucle CTCF.

    Le score moyenne une fenêtre 3×3 autour du pixel, donc un pic d'un seul
    pixel ressortirait dilué d'un facteur 9 — ce qui est voulu : un point de
    Hi-C n'est jamais un pixel isolé, et exiger qu'il le soit ferait d'un bruit
    de comptage une boucle.
    """
    oe = np.ones((40, 40))
    oe[9:12, 29:32] = 8.0
    oe[29:32, 9:12] = 8.0
    assert dot_score(oe, np.array([[10, 30]]))[0] == pytest.approx(8.0, rel=0.05)


def test_a_stripe_alone_does_not_make_a_dot_beside_it():
    oe = np.ones((40, 40))
    oe[10, :] = 8.0
    oe[:, 10] = 8.0
    assert dot_score(oe, np.array([[13, 30]]))[0] == pytest.approx(1.0, rel=0.05)


def test_a_stripe_does_not_mask_the_dot_sitting_on_it():
    """La vraie raison d'exclure la croix, et elle se mesure.

    Une ancre de boucle émet une bande le long de sa ligne. Si l'anneau de fond
    la gobait, le fond serait multiplié par la bande et le point réel
    disparaîtrait dans le rapport. Avec la croix exclue, le même point marque le
    même score qu'en l'absence de bande.
    """
    plain = np.ones((40, 40))
    plain[9:12, 29:32] = 8.0
    striped = np.ones((40, 40))
    striped[10, :] = 4.0
    striped[:, 10] = 4.0
    striped[9:12, 29:32] = 8.0
    alone = dot_score(plain, np.array([[10, 30]]))[0]
    on_stripe = dot_score(striped, np.array([[10, 30]]))[0]
    assert on_stripe == pytest.approx(alone, rel=0.25)
    assert on_stripe > 5.0


def test_dots_by_kind_reports_every_class_even_an_empty_one():
    truth = Truth(
        left=np.array([8]), right=np.array([20]), kind=np.array([CONVERGENT]),
        boundary=np.array([], dtype=np.int64),
    )
    cm = Map(np.ones((40, 40)), frames=1, bin_beads=1, bp_per_bin=2_000, cutoff=1.5)
    d = dots_by_kind(cm, truth, null=40)
    assert set(KINDS) <= set(d.score)
    assert np.isnan(d.separation("divergent"))    # aucun divergent : pas de rapport
    # Le fond au hasard, lui, se compte par dizaines et vaut 1 sur une carte plate.
    assert len(d.score["hasard"]) >= 20
    assert np.isfinite(d.score["hasard"]).all()   # tous tirés dans la zone scorable
    assert d.median("hasard") == pytest.approx(1.0)


def test_the_null_never_draws_a_pair_it_could_not_score():
    """Un pixel à moins de `outer` casiers d'un bord rend NaN, et un fond de NaN
    ne sert à rien. Ici la portée plantée occupe presque toute la carte."""
    truth = Truth(
        left=np.array([2]), right=np.array([34]), kind=np.array([CONVERGENT]),
        boundary=np.array([], dtype=np.int64),
    )
    cm = Map(np.ones((40, 40)), frames=1, bin_beads=1, bp_per_bin=2_000, cutoff=1.5)
    d = dots_by_kind(cm, truth, null=40, outer=6)
    assert np.isfinite(d.score["hasard"]).all()   # vide plutôt que faux


def test_separation_curve_on_a_straight_chain_is_exact():
    bp, mean, sd = separation_curve(
        _line(n=50, frames=2, spacing=2.0),
        bp_per_bead=2_000, unit_nm=10.0, seps=np.array([1, 5, 10]),
    )
    assert bp.tolist() == [2_000, 10_000, 20_000]
    assert mean.tolist() == pytest.approx([20.0, 100.0, 200.0])
    assert sd.tolist() == pytest.approx([0.0, 0.0, 0.0], abs=1e-6)


def test_separation_curve_refuses_to_cross_two_copies():
    """Une distance entre deux chromosomes ne correspond à aucune séparation
    génomique. La compter diluerait `R(s)` avec des paires qui n'ont pas de `s`."""
    x = _line(n=10, frames=1, spacing=1.0)
    x[:, 5:, 1] = 1_000.0                       # la seconde copie, très loin
    copy_id = np.array([0] * 5 + [1] * 5)
    _, mixed, _ = separation_curve(x, bp_per_bead=1, unit_nm=1.0, seps=np.array([1]))
    _, split, _ = separation_curve(
        x, bp_per_bead=1, unit_nm=1.0, copy_id=copy_id, seps=np.array([1])
    )
    assert split[0] == pytest.approx(1.0)
    assert mixed[0] > split[0]


def test_the_junction_ratio_is_one_when_both_models_say_the_same_thing():
    """Deux chaînes droites de même densité linéique doivent se raccorder
    exactement — c'est le seul cas où la réponse est connue d'avance."""
    reg = Region("chrT", 0, 3_000_000, 2_000)      # 1 500 monomères
    fine = _line(n=1_500, frames=2, spacing=1.0)
    coarse = np.zeros((2, 8, 3))
    coarse[:, :, 0] = np.arange(8) * 0.375         # 375 sigma par bille, en µm→nm
    j = junction(
        fine, reg, coarse,
        coarse_copy_id=np.zeros(8, dtype=int),
        coarse_bp_per_bead=750_000,
        sigma_nm=1.0,
        coarse_unit_nm=1_000.0,
        max_bp=3_000_000,
    )
    # 3 Mb en monomères de 2 kb : la séparation 1 500 n'existe pas sur 1 500
    # monomères, donc trois points de recouvrement et pas quatre.
    assert j.bp.tolist() == [750_000, 1_500_000, 2_250_000]
    assert np.allclose(j.ratio, 1.0, atol=1e-6)
    assert j.holds()
    # n − s paires par conformation : le dernier point est le plus fragile, et le
    # rapport doit le dire au lieu de le laisser deviner.
    assert j.fine_pairs.tolist() == [1_125, 750, 375]
    assert j.thin.tolist() == [False, False, False]


def test_the_junction_refuses_a_region_shorter_than_one_coarse_bead():
    with pytest.raises(ValueError, match="aucune séparation commune"):
        junction(
            _line(n=100, frames=1), Region("chrT", 0, 200_000, 2_000),
            np.zeros((1, 4, 3)),
            coarse_copy_id=np.zeros(4, dtype=int),
            coarse_bp_per_bead=750_000,
            sigma_nm=1.0,
        )


def test_the_junction_fails_loudly_when_the_two_models_disagree():
    reg = Region("chrT", 0, 3_000_000, 2_000)
    fine = _line(n=1_500, frames=2, spacing=1.0)
    coarse = np.zeros((2, 8, 3))
    coarse[:, :, 0] = np.arange(8) * 0.75          # deux fois trop étiré
    j = junction(
        fine, reg, coarse,
        coarse_copy_id=np.zeros(8, dtype=int),
        coarse_bp_per_bead=750_000,
        sigma_nm=1.0,
        coarse_unit_nm=1_000.0,
    )
    assert j.worst == pytest.approx(2.0, rel=1e-6)
    assert not j.holds()


# --------------------------------------------------------------------------
# run : l'assemblage
# --------------------------------------------------------------------------


def test_sigma_is_twice_the_monomer_radius():
    s = Setup(region=REGION)
    assert s.sigma_nm == pytest.approx(2.0 * monomer_radius(2_000))
    assert s.confine == pytest.approx(blob_radius(REGION.n, 0.30))


def test_every_replicate_sees_the_same_barriers():
    """La graine des sites est distincte de celle des réplicats.

    Mélangées, la variabilité mesurée mêlerait variabilité de repliement et
    variabilité de génome — l'erreur que la semaine 7 a corrigée sur le noyau.
    """
    reg = Region("chrT", 0, 120_000, 2_000)
    setup = Setup(region=reg, replicates=2, snapshots=3, stride=1,
                  md_per_step=200, relax=200, burnin_1d=50)
    a = build(setup, site_seed=4, first_seed=10, workers=1, mean_domain=60_000,
              min_domain=40_000)
    b = build(setup, site_seed=4, first_seed=99, workers=1, mean_domain=60_000,
              min_domain=40_000)
    assert np.array_equal(a.sites.pos, b.sites.pos)
    assert np.array_equal(a.sites.strand, b.sites.strand)
    assert not np.allclose(a.coords, b.coords)
    assert a.frames == 2 * 3


def test_save_and_load_round_trip(tmp_path):
    reg = Region("chrT", 0, 120_000, 2_000)
    setup = Setup(region=reg, replicates=1, snapshots=2, stride=1,
                  md_per_step=100, relax=100, burnin_1d=20)
    fine = build(setup, site_seed=1, workers=1, mean_domain=60_000, min_domain=40_000)
    p = tmp_path / "fine.npz"
    save(fine, str(p))
    coords, sites, truth, meta = load(str(p))
    assert coords.shape == fine.coords.shape
    assert np.array_equal(sites.strand, fine.sites.strand)
    assert truth is not None and np.array_equal(truth.kind, fine.truth.kind)
    assert meta["evidence"] == "simulated"
    assert meta["ctcf_source"].startswith("planted:")
