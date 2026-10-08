"""Les callers de conformation, validés contre une structure plantée.

C'est la seule étape du projet où « correct » a un sens. Sur des données réelles
on ne peut que constater un accord avec un autre outil ; ici on connaît la
réponse parce qu'on l'a écrite.

Les seuils sont volontairement stricts. Un caller qui passe de 100 % à 85 % de
rappel sur une structure plantée a un bug, pas un mauvais jour.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

pytest.importorskip("cooler", reason="pile Hi-C absente — voir docs/SETUP.md")
pytest.importorskip("cooltools", reason="pile Hi-C absente — voir docs/SETUP.md")

import cooler  # noqa: E402

from geno_pipeline.hic import (  # noqa: E402
    balance,
    compartments,
    gc_track,
    loops,
    match_pairs,
    match_positions,
    plant,
    tad_boundaries,
    write_cool,
)

logging.disable(logging.INFO)

N_BINS = 1_000
SEED = 3
WINDOW = 100_000


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    """Une matrice plantée, écrite en .cool réel et équilibrée. Construite une fois."""
    bins, pixels, truth = plant(n_bins=N_BINS, seed=SEED)
    path = write_cool(tmp_path_factory.mktemp("hic") / "planted.cool", bins, pixels)
    balance(cooler.Cooler(str(path)))
    return cooler.Cooler(str(path)), truth


# --------------------------------------------------------------------------
# L'instrument de mesure, d'abord
# --------------------------------------------------------------------------


def test_match_positions_counts_correctly():
    a = match_positions([10, 20, 30], [10, 20, 30], tol=0)
    assert (a.recall, a.precision, a.f1) == (1.0, 1.0, 1.0)

    b = match_positions([11, 20], [10, 20, 30], tol=1)
    assert b.n_matched == 2 and b.recall == pytest.approx(2 / 3) and b.precision == 1.0

    c = match_positions([500], [10, 20], tol=1)
    assert (c.recall, c.precision, c.f1) == (0.0, 0.0, 0.0)

    assert match_positions([], [10], tol=1).recall == 0.0
    assert match_positions([10], [], tol=1).n_true == 0


def test_match_pairs_requires_both_axes_close():
    true = np.array([[10, 50]])
    assert match_pairs(np.array([[10, 50]]), true, tol=0).recall == 1.0
    assert match_pairs(np.array([[11, 51]]), true, tol=1).recall == 1.0
    assert match_pairs(np.array([[11, 90]]), true, tol=1).recall == 0.0   # un seul axe proche


# --------------------------------------------------------------------------
# La structure plantée elle-même
# --------------------------------------------------------------------------


def test_planted_structure_is_coherent(planted):
    _clr, truth = planted

    sizes = np.diff(np.r_[0, truth.boundaries, truth.n_bins])
    assert sizes.min() >= 20, "un TAD de longueur nulle fabrique un positif inappelable"
    assert truth.boundaries.max() < truth.n_bins, "la fin du chromosome n'est pas une frontière"

    # Les transitions de compartiment doivent tomber sur des frontières de TAD.
    switches = np.flatnonzero(np.diff(truth.compartment)) + 1
    assert set(switches.tolist()) <= set(truth.boundaries.tolist())

    # Les boucles sont aux coins des TADs.
    for a, b in truth.loops:
        assert truth.domain[a] == truth.domain[b]
        assert a in truth.boundaries or a == 0


def test_gc_track_separates_the_compartments(planted):
    _clr, truth = planted
    assert truth.gc[truth.compartment == 1].mean() > truth.gc[truth.compartment == -1].mean()


# --------------------------------------------------------------------------
# Équilibrage
# --------------------------------------------------------------------------


def test_balancing_flattens_the_marginals(planted):
    clr, _truth = planted
    raw = clr.matrix(balance=False)[:]
    bal = clr.matrix(balance=True)[:]

    cv = lambda m: np.nanstd(np.nansum(m, axis=1)) / np.nanmean(np.nansum(m, axis=1))  # noqa: E731
    assert cv(bal) < cv(raw) / 5, "ICE doit aplatir les marginales, très nettement"
    assert cv(bal) < 0.05


def test_balancing_recovers_the_planted_bias(planted):
    clr, truth = planted
    w = clr.bins()["weight"][:].to_numpy()
    ok = np.isfinite(w)
    # Le poids corrige le biais : il doit varier à l'inverse de celui-ci.
    r = np.corrcoef(np.log(w[ok]), -np.log(truth.bias[ok]))[0, 1]
    assert r > 0.9, f"les poids ICE ne suivent pas le biais planté (r={r:.2f})"


# --------------------------------------------------------------------------
# Compartiments
# --------------------------------------------------------------------------


def test_eigenvector_recovers_the_partition(planted):
    clr, truth = planted
    e1 = compartments(clr).set_index("start")["E1"].to_numpy()
    ok = np.isfinite(e1)
    called = np.where(e1[ok] > 0, 1, -1)
    # Le signe est arbitraire sans phasage : on mesure la partition, pas l'étiquette.
    accuracy = max((called == truth.compartment[ok]).mean(), (called != truth.compartment[ok]).mean())
    assert accuracy > 0.95, f"partition mal retrouvée ({accuracy:.1%})"


def test_phasing_track_fixes_the_arbitrary_sign(planted):
    clr, truth = planted
    table = compartments(clr, phasing_track=gc_track(truth))
    e1 = table["E1"].to_numpy()
    ok = np.isfinite(e1)

    assert table.attrs["phased"] is True
    signed = np.where(e1[ok] > 0, 1, -1)
    assert (signed == truth.compartment[ok]).mean() > 0.95, "A doit sortir du bon côté"
    assert (table.loc[ok & (truth.compartment == 1), "compartment"] == "A").mean() > 0.95


# --------------------------------------------------------------------------
# TADs
# --------------------------------------------------------------------------


def _called_boundaries(clr, window):
    table = tad_boundaries(clr, window)
    return np.flatnonzero(table["is_boundary"].fillna(False).to_numpy())


def test_insulation_recovers_the_planted_boundaries(planted):
    clr, truth = planted
    a = match_positions(_called_boundaries(clr, WINDOW), truth.boundaries, tol=1)
    assert a.recall > 0.9, str(a)
    assert a.precision > 0.9, str(a)


def test_window_must_be_smaller_than_the_tad(planted):
    """Une fenêtre de la taille du TAD enjambe la frontière et lisse le minimum.

    C'est le réglage qu'on se serait trompé à copier d'un tutoriel sans regarder
    l'échelle de ses propres domaines.
    """
    clr, truth = planted
    tight = match_positions(_called_boundaries(clr, 100_000), truth.boundaries, tol=1)
    wide = match_positions(_called_boundaries(clr, 600_000), truth.boundaries, tol=1)
    assert tight.recall > wide.recall + 0.3, f"serrée {tight} / large {wide}"


def test_boundary_calls_land_within_one_bin(planted):
    """Elles ne tombent pas au bin exact, et c'est attendu : l'insulation est lissée.

    Toute comparaison de frontières, y compris celle à venir contre Rao 2014,
    doit donc se faire avec une tolérance — sinon on mesure du bruit.
    """
    clr, truth = planted
    called = _called_boundaries(clr, WINDOW)
    exact = match_positions(called, truth.boundaries, tol=0)
    loose = match_positions(called, truth.boundaries, tol=1)
    assert loose.recall > exact.recall


# --------------------------------------------------------------------------
# Boucles
# --------------------------------------------------------------------------


def test_loops_recover_the_planted_anchors(planted):
    clr, truth = planted
    table = loops(clr)
    pairs = np.c_[
        table["start1"].to_numpy() // clr.binsize,
        table["start2"].to_numpy() // clr.binsize,
    ]
    a = match_pairs(pairs, truth.loops, tol=2)
    assert a.recall > 0.8, str(a)
    assert a.precision > 0.9, str(a)
