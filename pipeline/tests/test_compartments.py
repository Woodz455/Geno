"""Semaine 10 — compartiments lus dans la carte de contacts d'un ensemble.

Même principe qu'en semaine 3 : on plante un damier connu dans un ensemble jouet, et on
vérifie que l'appel le retrouve ; puis on donne un ensemble **sans** damier, et on vérifie
que les trois mesures le disent. Un appelant qui trouverait des compartiments partout ne
passerait pas le second test.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("cooltools", reason="pile Hi-C absente — voir docs/SETUP.md")

from geno_pipeline.ensemble.compartments import compute, haploid_bins  # noqa: E402

BLOCK = 8          # billes par bloc de compartiment
PER_COPY = 96      # billes par copie


def toy(planted: bool, n_struct: int = 24, seed: int = 0):
    """Deux chromosomes diploïdes plus un X (Xa, Xi). Blocs A/B alternés ; les blocs B
    portent la piste LAD. Planté : A et B occupent deux régions distinctes de l'espace,
    donc A touche A et B touche B. Témoin : tout est mélangé dans une même sphère."""
    labels = ("chr1:a", "chr1:b", "chr2:a", "chr2:b", "chrX:Xa", "chrX:Xi")
    copy_id = np.repeat(np.arange(len(labels)), PER_COPY)
    # Comme le vrai génome : billes égales dans un chromosome, pas d'un chromosome à l'autre.
    # La première version du test avait des billes toutes égales et a laissé passer le refus
    # de cooltools (« mismatch between track and cooler bin size ») sur les vraies données.
    size = np.array([1000, 1000, 1013, 1013, 994, 994])[np.repeat(np.arange(len(labels)), PER_COPY)]
    start = np.tile(np.arange(PER_COPY), len(labels)) * size
    end = start + size
    is_b = (np.tile(np.arange(PER_COPY), len(labels)) // BLOCK) % 2 == 1
    lad = np.where(is_b, 0.9, 0.1)
    rng = np.random.default_rng(seed)
    n = len(copy_id)

    def ball(k, centre, r):
        u = rng.normal(size=(k, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        return centre + u * r * np.cbrt(rng.uniform(size=(k, 1)))

    coords = np.empty((n_struct, n, 3))
    for s in range(n_struct):
        if planted:
            coords[s, ~is_b] = ball(int((~is_b).sum()), np.array([-6.0, 0, 0]), 7.0)
            coords[s, is_b] = ball(int(is_b.sum()), np.array([6.0, 0, 0]), 7.0)
        else:
            coords[s] = ball(n, np.zeros(3), 9.0)
    radius = np.full(n, 0.6)
    return coords, radius, copy_id, labels, start, end, lad, is_b


def run(planted: bool):
    coords, radius, copy_id, labels, start, end, lad, is_b = toy(planted)
    return compute(coords, radius, copy_id, labels, start, end, lad, 15.0, cutoff=1.5), is_b, labels, copy_id


def test_a_planted_checkerboard_is_found_and_oriented_by_the_lad_track():
    c, is_b, labels, copy_id = run(planted=True)
    xa = np.array([labels[k] != "chrX:Xi" for k in copy_id])
    ok = np.isfinite(c.e1) & xa
    agree = ((c.e1[ok] > 0) == ~is_b[ok]).mean()
    assert agree > 0.95, agree
    assert c.strength > 2.0, c.strength
    assert c.split_half_r > 0.8, c.split_half_r
    assert c.r_lad < -0.8


def test_without_a_checkerboard_the_cross_validated_measures_say_so():
    """Le témoin : un E1 se calcule toujours, mais il ne veut rien dire, et ça se voit —
    à condition de croiser. Mesuré sur la carte dont il est tiré, le saddle d'un ensemble
    sans compartiment vaut ~2,8 : c'est ce biais que ce test garde en vue."""
    c, *_ = run(planted=False)
    assert c.strength < 1.3, c.strength
    assert c.strength_in_sample > 1.5, c.strength_in_sample
    assert abs(c.split_half_r) < 0.5, c.split_half_r


def test_the_inactive_x_gets_no_compartment():
    c, _, labels, copy_id = run(planted=True)
    xi = np.array([labels[k] == "chrX:Xi" for k in copy_id])
    xa = np.array([labels[k] == "chrX:Xa" for k in copy_id])
    assert np.isnan(c.e1[xi]).all()
    assert np.isfinite(c.e1[xa]).mean() > 0.9


def test_homologs_share_bins_and_a_different_cut_is_refused():
    labels = ("chr1:a", "chr1:b")
    copy_id = np.repeat([0, 1], 4)
    start = np.tile([0, 10, 20, 30], 2)
    end = start + 10
    bin_of, bins, xi = haploid_bins(labels, copy_id, start, end)
    assert bin_of.tolist() == [0, 1, 2, 3, 0, 1, 2, 3]
    assert len(bins) == 4 and not xi.any()
    start2 = start.copy()
    start2[5] = 11
    with pytest.raises(ValueError, match="découpée"):
        haploid_bins(labels, copy_id, start2, end)
