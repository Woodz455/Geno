"""Noyau diploïde complet : ce qui est imposé, ce qui est hérité, ce qui est mesuré.

Un test ici ne vaut que s'il distingue les trois. « Aucune interpénétration » est
une condition que le solveur impose : le tester, c'est tester le solveur. La
territorialité est héritée de l'initialisation : la tester, c'est vérifier que le
recuit ne la détruit pas, pas qu'elle est vraie. Et la stratification radiale
vient d'un terme du modèle nourri par une piste LAD synthétique : rien ici ne dit
quoi que ce soit sur un vrai noyau.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

pytest.importorskip("scipy", reason="pile de conformation absente — voir docs/SETUP.md")

from geno_pipeline.nucleus import (  # noqa: E402
    build,
    capacity_at,
    chrom_sizes,
    gm12878,
    periphery,
    save,
    segment,
    territoriality,
)
from geno_pipeline.nucleus.beads import PACKING_2D  # noqa: E402
from geno_pipeline.nucleus.karyotype import Copy, Karyotype  # noqa: E402
from geno_pipeline.nucleus.pack import (  # noqa: E402
    System,
    sealed,
    _radial_bias,
    _tighten,
    quality,
    relax,
)

BP_PER_BEAD = 750_000


@pytest.fixture(scope="module")
def karyotype():
    return gm12878()


@pytest.fixture(scope="module")
def beads(karyotype):
    return segment(karyotype, bp_per_bead=BP_PER_BEAD)


def tiny() -> Karyotype:
    """Six copies de 40 Mb — de quoi construire un vrai noyau miniature en 2 s."""
    copies = tuple(
        Copy(f"chr{i // 2 + 1}", i % 2, 40_000_000, f"chr{i // 2 + 1}:{'ab'[i % 2]}")
        for i in range(6)
    )
    return Karyotype("test 6 copies", copies, "test", "builtin")


@pytest.fixture(scope="module")
def small():
    return build(
        tiny(),
        bp_per_bead=2_000_000,
        nuclear_radius=1.0,
        steps=300,
        polish=1_200,
        territory_steps=200,
        seed=4,
    )


# --------------------------------------------------------------------------
# Caryotype : ce qu'on modélise, et d'où viennent ses longueurs
# --------------------------------------------------------------------------


def test_gm12878_is_46_XX(karyotype):
    assert karyotype.n_copies == 46
    chroms = [c.chrom for c in karyotype.copies]
    assert chroms.count("chrX") == 2, "GM12878 est une lignée féminine"
    assert "chrY" not in chroms
    for i in range(1, 23):
        assert chroms.count(f"chr{i}") == 2, f"chr{i} n'est pas en deux exemplaires"


def test_the_two_X_are_named_apart(karyotype):
    labels = [c.label for c in karyotype.copies if c.chrom == "chrX"]
    assert sorted(labels) == ["chrX:Xa", "chrX:Xi"]
    assert len(set(c.label for c in karyotype.copies)) == 46


def test_total_is_a_diploid_genome(karyotype):
    assert 5.9e9 < karyotype.total_bp < 6.2e9


def test_locked_file_wins_over_the_builtin_table(tmp_path):
    path = tmp_path / "hg38.chrom.sizes"
    path.write_text("chr1\t100\nchr1_KI270706v1_random\t50\nchrX\t70\n", encoding="utf-8")
    sizes, provenance = chrom_sizes(path)
    assert sizes == {"chr1": 100, "chrX": 70}, "les scaffolds ne sont pas l'assemblage primaire"
    assert provenance == str(path)


def test_builtin_provenance_is_announced(tmp_path):
    _, provenance = chrom_sizes(tmp_path / "jamais-récupéré")
    assert provenance == "builtin", "un modèle sur table interne doit le dire"


# --------------------------------------------------------------------------
# Billes : découpage, rayons, piste LAD
# --------------------------------------------------------------------------


def test_bead_count_lands_in_the_roadmap_band(beads):
    assert 6_000 <= beads.n <= 10_000, f"{beads.n} billes, hors de la fourchette semaine 6"


def test_bonds_never_cross_a_chromosome(beads):
    b = beads.bonds
    assert (beads.copy_id[b[:, 0]] == beads.copy_id[b[:, 1]]).all()
    # 46 chaînes séparées, donc 46 liaisons de moins qu'une chaîne unique.
    assert len(b) == beads.n - len(beads.labels)


def test_bead_volume_is_exactly_phi_of_the_nucleus(beads):
    ratio = (beads.radius**3).sum() / beads.nuclear_radius**3
    assert ratio == pytest.approx(beads.phi, rel=1e-9)


def test_radius_follows_the_cube_root_of_length(beads):
    density = (beads.end - beads.start) / beads.radius**3
    assert density.std() / density.mean() < 1e-9, "la densité locale doit être constante"


def test_lad_is_a_fraction_and_covers_a_plausible_share(beads):
    assert beads.lad.min() >= 0.0 and beads.lad.max() <= 1.0
    coverage = (beads.lad * beads.radius**3).sum() / (beads.radius**3).sum()
    assert 0.20 < coverage < 0.50, f"couverture LAD {coverage:.1%} hors du plausible"


def test_lad_is_blocky_not_drawn_bead_by_bead(beads):
    """Un LAD est un bloc. Une piste tirée bille par bille serait inutilisable :
    c'est l'autocorrélation qui décide si les billes LAD d'une copie peuvent
    atteindre la lamina ensemble."""
    inside = beads.copy_id[:-1] == beads.copy_id[1:]
    a, b = beads.lad[:-1][inside], beads.lad[1:][inside]
    assert float(np.corrcoef(a, b)[0, 1]) > 0.30


def test_bond_limit_leaves_room_for_contact(beads):
    b = beads.bonds
    contact = beads.radius[b[:, 0]] + beads.radius[b[:, 1]]
    assert (beads.bond_limits() > contact).all(), "une liaison ne doit pas exclure le contact"


def test_radial_targets_preserve_uniform_density(beads):
    """Les cibles suivent q^(1/3) : le profil visé est celui d'une sphère uniforme.

    Prendre la fraction LAD comme rayon visé tasserait tout le noyau vers son
    centre — le biais doit trier, pas comprimer.
    """
    t = beads.radial_targets()
    assert t.min() > 0.0 and t.max() < 1.0
    assert float((t**3).mean()) == pytest.approx(0.5, abs=0.01)


def test_radial_targets_are_ordered_by_lad(beads):
    t = beads.radial_targets()
    q = np.quantile(beads.lad, [0.25, 0.75])
    assert t[beads.lad >= q[1]].mean() > t[beads.lad <= q[0]].mean() + 0.2


def test_lamina_capacity_falls_with_resolution(karyotype):
    caps = [capacity_at(karyotype.total_bp, bp)[2] for bp in (3_000_000, 750_000, 100_000, 10_000)]
    assert caps == sorted(caps, reverse=True), "un modèle plus fin a moins de place à la lamina"
    assert caps[-1] < 0.15, "à 10 kb la borne doit être bien sous la couverture LAD"


def test_capacity_matches_the_closed_form(karyotype, beads):
    n, r_nm, cap = capacity_at(karyotype.total_bp, BP_PER_BEAD)
    r, R = r_nm / 1000.0, 5.0
    assert cap == pytest.approx(PACKING_2D * 4.0 * (R - r) ** 2 / r**2 / n, rel=1e-6)
    assert r_nm == pytest.approx(beads.radius.mean() * 1000, rel=0.01)


# --------------------------------------------------------------------------
# Solveur : les conditions qu'il impose
# --------------------------------------------------------------------------


def test_relax_separates_a_pile_of_spheres():
    """`tol` est passé explicitement : le polissage s'arrête dès qu'il est tenu,
    donc un test qui exigerait mieux que la tolérance sous laquelle il tourne
    mesurerait le réglage par défaut et non le solveur."""
    rng = np.random.default_rng(0)
    n = 300
    sys = System(x=rng.normal(0, 0.3, (n, 3)), r=np.full(n, 0.12), R=3.0)
    q, _ = relax(sys, steps=200, polish=1_500, inflate_from=0.4, tol=0.002, seed=1)
    assert q.max_overlap < 0.005, str(q)
    assert q.outside == 0


def test_confinement_is_hard():
    rng = np.random.default_rng(1)
    n = 400
    sys = System(x=rng.normal(0, 2.0, (n, 3)), r=np.full(n, 0.1), R=1.5)
    q, _ = relax(sys, steps=200, polish=800, inflate_from=0.4, seed=2)
    assert q.outside == 0
    assert (np.linalg.norm(sys.x, axis=1) + sys.r <= 1.5 + 1e-9).all()


def test_quality_sees_a_planted_overlap():
    sys = System(
        x=np.array([[0.0, 0, 0], [0.5, 0, 0]]), r=np.array([0.5, 0.5]), R=10.0
    )
    assert quality(sys).max_overlap == pytest.approx(0.5)
    assert quality(sys).n_overlapping == 1


def test_bond_is_one_sided():
    """Une liaison limite l'allongement et rien d'autre : comprimée, elle ne tire pas.

    Une contrainte bilatérale se bat contre le volume exclu dans les replis
    serrés, et c'est le volume exclu qui cède — or c'est lui la condition.
    """
    bonds, limit = np.array([[0, 1]]), np.array([1.0])

    short = System(x=np.array([[0.0, 0, 0], [0.4, 0, 0]]), r=np.full(2, 0.2),
                   R=10.0, bonds=bonds, bond_len=limit)
    before = short.x.copy()
    _tighten(short, 0.5)
    assert np.allclose(short.x, before), "une liaison comprimée ne doit pas tirer"

    long = System(x=np.array([[0.0, 0, 0], [2.0, 0, 0]]), r=np.full(2, 0.2),
                  R=10.0, bonds=bonds, bond_len=limit)
    _tighten(long, 1.0)
    assert np.linalg.norm(long.x[1] - long.x[0]) == pytest.approx(1.0)


def test_radial_bias_converges_instead_of_accumulating():
    """Régression. Un biais radial doit être un **rappel**, pas une poussée.

    Une poussée s'accumule sur tous les pas du recuit : son effet dépend alors du
    nombre de pas et non du modèle. Mesuré à l'époque, 500 pas d'une poussée même
    faible plaquaient toutes les billes LAD contre l'enveloppe et y créaient une
    croûte bloquée — chevauchement résiduel 21 %, contre 5 % sans poussée du tout.
    Ce test échoue pour toute formulation qui accumule.
    """

    def run(steps: int) -> float:
        sys = System(x=np.array([[1.0, 0.0, 0.0]]), r=np.array([0.1]), R=10.0,
                     outward=np.array([0.5]))
        for _ in range(steps):
            _radial_bias(sys, 0.1)
        return float(np.linalg.norm(sys.x[0]))

    target = (10.0 - 0.1) * 0.5
    assert run(300) == pytest.approx(target, abs=1e-6)
    assert run(3_000) == pytest.approx(target, abs=1e-6), "dix fois plus de pas, même réponse"


# --------------------------------------------------------------------------
# Mesures
# --------------------------------------------------------------------------


def test_an_octahedral_cage_is_a_fixed_point_that_only_a_shake_breaks():
    """Régression. Six billes en contact mutuel forment une cage rigide.

    Trois liaisons tombées sur les trois diagonales d'un octaèdre ne peuvent plus
    se raccourcir : raccourcir une diagonale écarterait les quatre billes de
    l'équateur, ce que les deux autres liaisons interdisent. Le rapport
    diagonale/arête valant √2, la tension se fige à `√2 · 0,98 / 1,15 − 1`, soit
    **+20,4 %** — la constante exacte relevée sur cinq configurations et trois
    graines du vrai noyau.

    Ce n'est pas un puits peu profond : on n'en sort pas par agitation, il faut
    casser la cage.
    """
    r, tight = 0.5, 0.98
    half = tight * 2.0 * r / np.sqrt(2.0)          # demi-diagonale de l'octaèdre
    x = np.array([[half, 0, 0], [-half, 0, 0], [0, half, 0],
                  [0, -half, 0], [0, 0, half], [0, 0, -half]])
    bonds = np.array([[0, 1], [2, 3], [4, 5]])     # une liaison par diagonale
    limit = np.full(3, 1.15 * 2.0 * r)

    def cage() -> System:
        return System(x=x.copy(), r=np.full(6, r), R=10.0, bonds=bonds, bond_len=limit)

    start = quality(cage())
    assert start.bond_stretch == pytest.approx(np.sqrt(2) * tight / 1.15 - 1.0, abs=1e-6)
    assert start.max_overlap == pytest.approx(1.0 - tight, abs=1e-6)

    stuck, _ = relax(cage(), steps=0, polish=800, polish_t=0.0, shake_rounds=0, seed=0)
    assert stuck.bond_stretch > 0.15, f"la cage aurait dû tenir : {stuck}"

    broken, _ = relax(cage(), steps=0, polish=800, polish_t=0.0, shake_rounds=6, seed=0)
    assert broken.acceptable(0.01), f"la secousse n'a pas cassé la cage : {broken}"


def test_the_sealing_condition_is_s_above_half_the_stretch():
    """Pendant l'inflation, le tube de la chaîne doit rester étanche.

    À l'échelle `s`, deux billes liées s'excluent au rayon `s·r` mais restent
    séparées d'au plus `stretch·(r_i + r_j)` : une troisième bille passe entre
    elles dès que `s ≤ stretch / 2`. Mesuré à s = 0,50 avec stretch = 1,15 : une
    liaison bloquée 20 % au-delà de sa limite et 5 % de chevauchement résiduel là
    où tout le reste était à 0,04 %.
    """
    assert not sealed(0.50, 1.15)
    assert not sealed(0.575, 1.15)
    assert sealed(0.58, 1.15)
    assert sealed(0.65, 1.15)


def test_the_default_build_stays_sealed(small):
    assert small.sealed, "les valeurs par défaut doivent laisser le tube étanche"


def test_territoriality_is_one_when_everything_is_mixed():
    rng = np.random.default_rng(3)
    x = rng.normal(0, 1, (1_200, 3))
    copy_id = rng.integers(0, 6, 1_200)
    t = territoriality(x, copy_id)
    assert 0.85 < t.index < 1.15, f"indice {t.index:.2f} sur un mélange parfait"


def test_territoriality_is_large_when_chromosomes_are_apart():
    rng = np.random.default_rng(4)
    centres = np.array([[0, 0, 0], [8, 0, 0], [0, 8, 0], [0, 0, 8], [8, 8, 0], [8, 0, 8]])
    copy_id = np.repeat(np.arange(6), 200)
    x = centres[copy_id] + rng.normal(0, 1, (1_200, 3))
    assert territoriality(x, copy_id).index > 5.0


def test_uniform_points_give_an_enrichment_of_one():
    """La part « densité uniforme » doit être exactement la part de volume de la
    coquille — c'est elle qui rend l'enrichissement mesuré interprétable."""
    rng = np.random.default_rng(5)
    n, R, r = 60_000, 5.0, 0.15
    u = rng.random(n) ** (1 / 3) * (R - r)
    d = rng.normal(0, 1, (n, 3))
    x = d / np.linalg.norm(d, axis=1)[:, None] * u[:, None]
    p = periphery(x, np.full(n, r), rng.random(n), R, capacity=1_000)
    assert p.enrichment == pytest.approx(1.0, abs=0.05), str(p)
    assert p.uniform_share == pytest.approx(1.5 * r / (R - r), rel=0.05)


# --------------------------------------------------------------------------
# Noyau complet
# --------------------------------------------------------------------------


def test_small_nucleus_has_no_interpenetration(small):
    assert small.final.acceptable(0.01), str(small.final)
    assert small.final.outside == 0


def test_chains_stay_inside_the_nucleus(small):
    b = small.beads
    assert (np.linalg.norm(small.x, axis=1) + b.radius <= b.nuclear_radius + 1e-9).all()


def test_chains_never_break(small):
    b = small.beads
    bonds = b.bonds
    d = np.linalg.norm(small.x[bonds[:, 0]] - small.x[bonds[:, 1]], axis=1)
    assert (d <= b.bond_limits() * 1.05).all(), "une liaison a cassé"


def test_anneal_preserves_the_territories_it_inherited(small):
    """La territorialité est une **entrée** — une relaxation ne fait jamais se
    croiser deux chaînes. Le seul énoncé vérifiable est que le recuit la garde."""
    assert small.before.index > 3.0, "l'initialisation n'a pas fait de territoires"
    assert small.after.index > 0.5 * small.before.index, (
        f"le recuit a dissous les territoires : {small.before.index:.1f} → "
        f"{small.after.index:.1f}"
    )


def test_saved_nucleus_carries_its_provenance(small, tmp_path):
    sidecar = save(small, tmp_path / "noyau.npz")
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    assert meta["evidence"] == "simulated"
    assert meta["lad_source"] == "synthetic"
    assert meta["chrom_sizes_provenance"] == "builtin"
    assert "simulée" in meta["warning"]

    data = np.load(tmp_path / "noyau.npz")
    assert data["x"].shape == (small.n, 3)
    assert set(data.files) >= {"x", "radius", "copy_id", "start", "end", "lad"}
