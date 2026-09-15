"""Ensembles : un génome, N repliements, et les mesures qui ont un sens sur eux.

La ligne de partage de la semaine est géométrique. Deux noyaux recuits séparément
ne partagent aucun repère, donc une variance par bille en x, y, z serait un
nombre sans objet. Les tests d'ici vérifient d'abord ça — que le module mesure
l'absence de repère au lieu de la supposer — puis que chaque statistique
invariante par rotation rend la bonne réponse sur une entrée dont on connaît la
réponse.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy", reason="pile de conformation absente — voir docs/SETUP.md")
pytest.importorskip("zarr", reason="zarr absent — voir docs/SETUP.md")

from geno_pipeline.ensemble import (  # noqa: E402
    contact_pairs,
    contacts,
    create,
    damid,
    frame,
    generate,
    homolog_map,
    medoid,
    pending,
    put,
    read,
    read_bedgraph,
    radial_by_quartile,
    reproducibility,
    self_consistency,
    slope,
)
from geno_pipeline.nucleus import gm12878, segment  # noqa: E402
from geno_pipeline.nucleus.pack import Quality  # noqa: E402


@pytest.fixture(scope="module")
def beads():
    return segment(gm12878(), bp_per_bead=3_000_000)


def line(n: int, radius: float, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    """Chaîne rectiligne : on sait exactement quelles paires sont en contact."""
    x = np.zeros((n, 3))
    x[:, 0] = np.arange(n) * spacing
    return x, np.full(n, radius)


# --------------------------------------------------------------------------
# Magasin
# --------------------------------------------------------------------------


def test_store_round_trip(beads, tmp_path):
    path = tmp_path / "e.zarr"
    root = create(path, beads, 3, {"karyotype": "test", "lad_seed": 0})
    rng = np.random.default_rng(0)
    written = {}
    for k in range(3):
        coords = rng.normal(0, 1, (beads.n, 3))
        written[k] = coords
        put(root, k, 100 + k, coords, Quality(0.004, 0.001, 7, 0, 0.002, shakes=k))

    e = read(path)
    assert e.n_structures == 3 and e.n_beads == beads.n
    assert list(e.seeds) == [100, 101, 102]
    assert np.allclose(e.coords[1], written[1], atol=1e-6)
    assert np.array_equal(e.beads.copy_id, beads.copy_id)
    assert np.allclose(e.beads.lad, beads.lad, atol=1e-6)
    assert e.meta["evidence"] == "simulated"
    assert e.meta["lad_source"] == "synthetic"
    assert "simulé" in e.meta["warning"]


def test_unfinished_structures_are_not_returned(beads, tmp_path):
    """Une génération interrompue ne doit pas rendre des zéros pour des structures
    jamais produites — un zéro est une position, pas une absence."""
    path = tmp_path / "e.zarr"
    root = create(path, beads, 5, {})
    put(root, 0, 7, np.ones((beads.n, 3)), Quality(0.001, 0.0, 0, 0, 0.0))
    put(root, 3, 9, np.ones((beads.n, 3)) * 2, Quality(0.001, 0.0, 0, 0, 0.0))

    e = read(path)
    assert e.n_structures == 2
    assert list(e.seeds) == [7, 9]


# --------------------------------------------------------------------------
# Le génome doit être le même partout
# --------------------------------------------------------------------------


def test_conformation_seed_does_not_move_the_genome():
    """Régression sur la raison d'être de `lad_seed`.

    Faire varier la seule graine ferait varier la piste LAD d'une structure à
    l'autre : on mesurerait alors la variabilité de deux cents génomes différents
    en croyant mesurer celle d'un repliement.
    """
    a = segment(gm12878(), bp_per_bead=3_000_000, seed=1)
    b = segment(gm12878(), bp_per_bead=3_000_000, seed=2)
    assert not np.allclose(a.lad, b.lad), "la graine doit bien changer le génome"

    same = segment(gm12878(), bp_per_bead=3_000_000, seed=1)
    assert np.array_equal(a.lad, same.lad), "à graine égale, génome égal"


# --------------------------------------------------------------------------
# Il n'y a pas de repère commun
# --------------------------------------------------------------------------


def test_frame_finds_a_shared_frame_when_there_is_one():
    """Témoin positif : des copies tournées de la même structure **partagent** un repère."""
    rng = np.random.default_rng(0)
    base = rng.normal(0, 1, (400, 3))
    stack = []
    for _ in range(6):
        q, _ = np.linalg.qr(rng.normal(0, 1, (3, 3)))
        stack.append(base @ q)
    f = frame(np.array(stack), nuclear_radius=1.0, pairs=6, seed=1)
    assert f.aligned < 1e-6, str(f)
    assert f.kept < 0.01


def test_frame_finds_none_between_independent_structures():
    """Témoin négatif : des nuages indépendants n'en partagent aucun, et l'alignement
    optimal ne rattrape presque rien."""
    rng = np.random.default_rng(1)
    stack = rng.normal(0, 1, (6, 400, 3))
    f = frame(stack, nuclear_radius=1.0, pairs=6, seed=2)
    assert f.kept > 0.9, str(f)


# --------------------------------------------------------------------------
# Reproductibilité
# --------------------------------------------------------------------------


def test_icc_is_one_when_every_structure_agrees():
    rng = np.random.default_rng(2)
    profile = rng.random(500)
    radial = np.tile(profile, (30, 1))
    rep = reproducibility(radial, 5.0, pairs=20)
    assert rep.icc > 0.999 and rep.pairwise_r > 0.999
    assert rep.within_nm < 1e-6


def test_icc_is_zero_when_depth_is_pure_draw():
    rng = np.random.default_rng(3)
    radial = rng.random((30, 500))
    rep = reproducibility(radial, 5.0, pairs=20)
    assert abs(rep.icc) < 0.05, str(rep)
    assert abs(rep.pairwise_r) < 0.10


def test_medoid_picks_the_central_structure():
    rng = np.random.default_rng(4)
    centre = rng.random(300)
    radial = np.vstack([centre] + [centre + rng.normal(0, 0.3, 300) for _ in range(9)])
    index, total = medoid(radial)
    assert index == 0, f"médoïde {index}, distances {total.round(1)}"


def test_self_consistency_follows_the_lad_track():
    rng = np.random.default_rng(5)
    lad = rng.random(400)
    radial = np.tile(lad, (10, 1)) + rng.normal(0, 0.01, (10, 400))
    assert self_consistency(radial, lad) > 0.99


# --------------------------------------------------------------------------
# Contacts
# --------------------------------------------------------------------------


def test_homolog_map_reads_the_labels(beads):
    partner = homolog_map(beads.labels)
    for k, label in enumerate(beads.labels):
        chrom = label.split(":")[0]
        assert beads.labels[partner[k]].split(":")[0] == chrom
        assert partner[k] != k, f"{label} serait son propre homologue"


def test_homolog_map_leaves_an_unpaired_chromosome_alone():
    """Une trisomie casse l'appariement sans casser l'ordre : `index ^ 1` mentirait."""
    labels = ("chr21:a", "chr21:b", "chr21:c", "chrX:Xa", "chrX:Xi")
    partner = homolog_map(labels)
    assert list(partner[:3]) == [0, 1, 2], "trois copies n'ont pas d'homologue unique"
    assert list(partner[3:]) == [4, 3]


def test_contact_pairs_respects_the_cutoff():
    x, r = line(4, radius=0.5, spacing=1.0)      # contact exact à d = 1,0
    assert len(contact_pairs(x, r, cutoff=0.9)) == 0
    assert len(contact_pairs(x, r, cutoff=1.5)) == 3        # les 3 paires voisines
    assert len(contact_pairs(x, r, cutoff=2.5)) == 5        # + les 2 paires à 2 crans


def test_ps_normalisation_gives_exactly_one_where_everything_touches():
    """Sur une chaîne droite dont on connaît la portée, P(s) doit valoir 1 jusqu'à
    la portée et 0 après. C'est la normalisation par le nombre de paires
    *possibles* qui est testée : sans elle, P(s) décroîtrait toute seule parce
    qu'il y a moins de paires disponibles aux grandes séparations."""
    n = 60
    x, r = line(n, radius=0.5, spacing=1.0)
    coords = np.stack([x, x])                     # deux structures identiques
    copy_id = np.zeros(n, dtype=np.int64)
    c = contacts(
        coords, r, copy_id, ("chr1:a",), bp_per_bead=1e6, cutoff=2.5,
        fit_from=1e9, fit_to=2e9,                 # fenêtre vide : la pente reste indéfinie
    )
    p = dict(zip((c.separation / 1e6).round().astype(int), c.probability.round(9)))
    assert p[1] == 1.0 and p[2] == 1.0
    assert 3 not in p, "rien ne doit toucher au-delà de la portée"
    assert c.cis == 2 * ((n - 1) + (n - 2)), "cis est un total d'ensemble, pas un par-structure"
    assert c.per_structure == (n - 1) + (n - 2)
    assert c.trans == 0
    assert c.homolog == 0, "une copie sans homologue ne contacte pas d'homologue"
    assert np.isnan(c.exponent), "une fenêtre sans point ne doit pas produire de pente"


def test_trans_and_homolog_contacts_are_counted_apart():
    # Deux copies d'un même chromosome, superposées : tout contact est un homologue.
    n = 6
    x, r = line(n, radius=0.5, spacing=1.0)
    both = np.vstack([x, x + np.array([0.0, 0.2, 0.0])])
    coords = both[None, :, :]
    copy_id = np.repeat([0, 1], n)
    c = contacts(
        coords, np.concatenate([r, r]), copy_id, ("chr9:a", "chr9:b"),
        bp_per_bead=1e6, cutoff=1.5,
    )
    assert c.trans > 0
    assert c.homolog == c.trans, "les deux copies sont homologues, donc tous les trans le sont"


# --------------------------------------------------------------------------
# Le livrable qui attend son fichier
# --------------------------------------------------------------------------


def test_damid_recovers_a_planted_correlation():
    """La corrélation contre un DamID *mesuré* est le livrable de la semaine 7. Le
    réseau étant fermé, on la teste sur une piste construite exprès."""
    n = 40
    starts = np.arange(n) * 1_000_000
    ends = starts + 1_000_000
    copy_id = np.zeros(n, dtype=np.int64)
    radial = np.tile(np.linspace(0.1, 0.9, n), (5, 1))

    track = {"chr1": [(int(s), int(e), float(v))
                      for s, e, v in zip(starts, ends, np.linspace(0.1, 0.9, n))]}
    d = damid(radial, starts, ends, copy_id, ("chr1:a",), track, "planté")
    assert d.covered == n
    assert d.correlation > 0.999, str(d)


def test_damid_excludes_beads_the_track_does_not_cover():
    """Une bille non couverte est une absence, pas un zéro : la mettre à zéro
    fabriquerait de la corrélation là où il n'y a pas de mesure."""
    n = 10
    starts = np.arange(n) * 1_000_000
    ends = starts + 1_000_000
    copy_id = np.zeros(n, dtype=np.int64)
    radial = np.tile(np.linspace(0.1, 0.9, n), (3, 1))

    half = {"chr1": [(int(starts[k]), int(ends[k]), float(k)) for k in range(4)]}
    d = damid(radial, starts, ends, copy_id, ("chr1:a",), half, "partiel")
    assert d.covered == 4 and d.total == n

    d2 = damid(radial, starts, ends, copy_id, ("chr7:a",), half, "mauvais chromosome")
    assert d2.covered == 0 and np.isnan(d2.correlation)


def test_bedgraph_reader_skips_headers_and_refuses_junk(tmp_path):
    path = tmp_path / "lads.bedGraph"
    path.write_text(
        "track type=bedGraph name=LAD\n"
        "# un commentaire\n"
        "\n"
        "chr1\t0\t1000\t0.5\n"
        "chr1\t1000\t2000\t-1.25\n"
        "chr2\t0\t500\t2\n",
        encoding="utf-8",
    )
    track, rows = read_bedgraph(path)
    assert rows == 3
    assert track["chr1"] == [(0, 1000, 0.5), (1000, 2000, -1.25)]
    assert track["chr2"] == [(0, 500, 2.0)]

    bad = tmp_path / "bad.bedGraph"
    bad.write_text("chr1\t0\t1000\tNaN-ish\n", encoding="utf-8")
    with pytest.raises(ValueError, match="bad.bedGraph:1"):
        read_bedgraph(bad)

    short = tmp_path / "short.bedGraph"
    short.write_text("chr1\t0\t1000\n", encoding="utf-8")
    with pytest.raises(ValueError, match="quatre colonnes"):
        read_bedgraph(short)


def test_damid_runs_end_to_end_from_a_file(tmp_path):
    """Le chemin complet — fichier, lecture, corrélation — tourne. Il n'a jamais
    tourné sur un vrai DamID, faute de réseau, mais il ne reste rien à écrire."""
    n = 30
    starts = np.arange(n) * 500_000
    ends = starts + 500_000
    values = np.linspace(-2.0, 2.0, n)

    path = tmp_path / "damid.bedGraph"
    path.write_text(
        "track type=bedGraph\n"
        + "".join(f"chr1\t{s}\t{e}\t{v:.4f}\n" for s, e, v in zip(starts, ends, values)),
        encoding="utf-8",
    )
    track, rows = read_bedgraph(path)
    assert rows == n

    radial = np.tile(values, (4, 1)) * 0.1 + 0.5
    d = damid(radial, starts, ends, np.zeros(n, dtype=np.int64), ("chr1:a",), track, str(path))
    assert d.covered == n and d.correlation > 0.999


def test_generate_completes_an_existing_store_instead_of_wiping_it(beads, tmp_path):
    """Une série de deux cents structures coûte une demi-heure. Relancer la commande
    ne doit pas l'effacer — et `zarr.open_group(mode="w")` le ferait sans rien
    demander. Reprendre est donc le défaut, recréer se demande explicitement."""
    path = tmp_path / "e.zarr"
    root = create(path, beads, 3, {"karyotype": "test"})
    kept = np.full((beads.n, 3), 0.25)
    put(root, 0, 42, kept, Quality(0.002, 0.0, 0, 0, 0.001))

    # Magasin incomplet mais non vide : rien n'est perdu, et il reste 2 à produire.
    assert list(pending(path)) == [1, 2]

    # Même appel, aucune structure à produire au-delà : les données restent en place.
    generate(path, n_structures=3, workers=1, bp_per_bead=3_000_000)
    e = read(path)
    assert e.n_structures >= 1
    assert np.allclose(e.coords[0], kept, atol=1e-6)
    assert int(e.seeds[0]) == 42, "la structure déjà en magasin a été écrasée"


def test_generate_refuses_a_store_of_a_different_size(beads, tmp_path):
    path = tmp_path / "e.zarr"
    create(path, beads, 3, {})
    with pytest.raises(ValueError, match="--fresh"):
        generate(path, n_structures=7, workers=1)


def test_slope_recovers_a_planted_power_law():
    s = np.logspace(6, 8, 40)
    p = 3.0 * s ** (-1.25)
    assert slope(s, p, 1e6, 1e8) == pytest.approx(-1.25, abs=1e-9)
    assert np.isnan(slope(s, p, 1e9, 2e9)), "une fenêtre vide ne produit pas de pente"


def test_quartiles_survive_a_track_that_is_mostly_zero(beads):
    """Régression. Une bille de 750 kb sans le moindre LAD est fréquente, si bien
    que les quartiles *de valeur* valent tous deux zéro et que Q1 ressort vide —
    le rapport affichait `nan`. Le découpage par rang est toujours défini."""
    lad = np.zeros(400)
    lad[300:] = np.linspace(0.1, 1.0, 100)
    radial = np.tile(np.linspace(0.2, 0.9, 400), (5, 1))

    groups = radial_by_quartile(radial, lad)
    assert len(groups) == 4
    for mean, sd, lo, hi in groups:
        assert np.isfinite(mean) and np.isfinite(sd)
        assert lo <= hi
    assert groups[-1][0] > groups[0][0], "le quartile le plus LAD doit être le plus externe"


def test_contacts_reports_every_named_regime():
    n = 300
    x, r = line(n, radius=0.5, spacing=1.0)
    coords = x[None, :, :]
    c = contacts(
        coords, r, np.zeros(n, dtype=np.int64), ("chr1:a",), bp_per_bead=1e6, cutoff=3.5
    )
    assert [name for name, _ in c.regimes] == ["chaîne", "polymère", "territoire"]
    assert c.exponent == pytest.approx(
        slope(c.separation, c.probability, *c.fit_range), nan_ok=True
    )
