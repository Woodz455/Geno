"""Semaine 9 — le format `.g3d`, écrit en Python et relu partout.

Tout est construit sur des niveaux synthétiques : le format ne doit dépendre d'aucune
donnée générée par une autre semaine pour être vérifiable.
"""

from __future__ import annotations

import json
import struct

import numpy as np
import pytest

from geno_pipeline.export import CorruptError, Level, read, read_header, write
from geno_pipeline.export.build import coarsen, pose
from geno_pipeline.export.g3d import (
    CHAINS,
    MAGIC,
    PREAMBLE,
    decode,
    encode,
    octree,
    precision,
)


def chain_level(n_per_copy=(300, 200, 250), name="noyau", seed=0, spread=4000.0) -> Level:
    """Des copies en marche aléatoire de pas ~400 nm, des billes de 750 kb."""
    rng = np.random.default_rng(seed)
    pos, copy, start, end = [], [], [], []
    for c, k in enumerate(n_per_copy):
        walk = np.cumsum(rng.normal(0, 230, (k, 3)), axis=0) + rng.normal(0, spread, 3)
        pos.append(walk)
        copy.append(np.full(k, c))
        s = np.arange(k) * 750_000
        start.append(s)
        end.append(s + 750_000)
    n = sum(n_per_copy)
    return Level(
        name=name,
        bp_per_bead=750_000,
        positions_nm=np.concatenate(pos),
        copy=np.concatenate(copy),
        start=np.concatenate(start),
        end=np.concatenate(end),
        variability_nm=rng.uniform(100, 400, n),
        copies=[f"chr{c + 1}:a" for c in range(len(n_per_copy))],
        radius_nm=[166.8] * len(n_per_copy),
    )


def by_locus(copy, start):
    return np.lexsort((start, copy))


# --------------------------------------------------------------------------
# Filtres
# --------------------------------------------------------------------------


@pytest.mark.parametrize("chain", CHAINS)
@pytest.mark.parametrize("dtype", [np.uint8, np.uint16, np.uint32])
def test_every_filter_chain_is_lossless(chain, dtype):
    rng = np.random.default_rng(1)
    a = rng.integers(0, np.iinfo(dtype).max, 3 * 501, dtype=dtype)
    back = decode(encode(a, 3, chain), np.dtype(dtype).str.lstrip("<|"), 3, list(chain))
    assert np.array_equal(back, a)


def test_delta_wraps_around_instead_of_going_negative():
    """En u16, 3 − 5 vaut 65 534, et la somme cumulée retombe sur 3. Une suite qui
    descend ne doit donc ni lever, ni se corrompre."""
    a = np.array([5, 3, 65535, 0, 1], dtype=np.uint16)
    back = decode(encode(a, 1, ("delta", "deflate")), "u2", 1, ["delta", "deflate"])
    assert back.tolist() == a.tolist()


def test_delta_restarts_on_each_plane():
    """x, y et z sont trois plans : le delta ne doit pas enjamber la frontière entre eux."""
    a = np.array([10, 11, 12, 900, 901, 902, 5, 6, 7], dtype=np.uint16)
    assert decode(encode(a, 3, ("delta", "deflate")), "u2", 3, ["delta", "deflate"]).tolist() \
        == a.tolist()


# --------------------------------------------------------------------------
# Écriture et relecture
# --------------------------------------------------------------------------


def test_round_trip_is_exact_for_integers_and_bounded_for_positions(tmp_path):
    lv = chain_level()
    rep = write(tmp_path / "a.g3d", [lv])
    got = read(rep["path"]).levels["noyau"]

    o1, o2 = by_locus(lv.copy, lv.start), by_locus(got["copy"], got["start"])
    assert np.array_equal(got["copy"][o2], lv.copy[o1])
    assert np.array_equal(got["start"][o2], lv.start[o1])
    assert np.array_equal(got["end"][o2], lv.end[o1])

    err = np.abs(got["positions_nm"][o2] - lv.positions_nm[o1]).max()
    announced = got["meta"]["quant"]["max_error_nm"]
    assert err <= announced * (1 + 1e-9)


def test_the_announced_error_follows_the_uncertainty_rule(tmp_path):
    """1 % de l'incertitude médiane du niveau, sauf si la boîte l'interdit en 16 bits."""
    lv = chain_level()
    rep = write(tmp_path / "a.g3d", [lv], fraction=0.01)
    meta = read(rep["path"]).levels["noyau"]["meta"]
    assert meta["quant"]["max_error_nm"] == pytest.approx(0.01 * np.median(lv.variability_nm))


def test_a_coarser_rule_makes_a_smaller_file(tmp_path):
    lv = chain_level()
    fine = write(tmp_path / "f.g3d", [lv], fraction=0.0001)["bytes"]
    coarse = write(tmp_path / "c.g3d", [lv], fraction=0.05)["bytes"]
    assert coarse < fine


def test_precision_refuses_a_level_without_uncertainty():
    with pytest.raises(ValueError, match="incertitude"):
        precision(np.zeros(10), 0.01)


def test_variability_survives_within_the_same_error(tmp_path):
    lv = chain_level()
    rep = write(tmp_path / "a.g3d", [lv])
    got = read(rep["path"]).levels["noyau"]
    o1, o2 = by_locus(lv.copy, lv.start), by_locus(got["copy"], got["start"])
    err = np.abs(got["variability_nm"][o2] - lv.variability_nm[o1]).max()
    assert err <= got["meta"]["variability"]["step"] / 2 * (1 + 1e-9)


# --------------------------------------------------------------------------
# Disposition : tout préfixe est utile
# --------------------------------------------------------------------------


def test_first_frame_columns_come_first_and_contiguous(tmp_path):
    """Les colonnes du premier rendu de **tous** les chunks du niveau désigné tiennent
    avant `first_frame_end`, et rien d'autre ne s'y glisse."""
    big = chain_level(n_per_copy=(3000, 2500), name="noyau")
    small = chain_level(n_per_copy=(400,), name="apercu", seed=3)
    rep = write(tmp_path / "a.g3d", [small, big], first="noyau", max_beads=1000)
    blob = (tmp_path / "a.g3d").read_bytes()
    hdr, base, first_end = read_header(blob)
    assert first_end == rep["first_frame_end"]

    level = next(lv for lv in hdr["levels"] if lv["name"] == "noyau")
    assert len(level["chunks"]) > 1
    ends = []
    for c in level["chunks"]:
        for name in ("copy", "position"):
            spec = c["columns"][name]
            ends.append(base + spec["offset"] + spec["length"])
    assert max(ends) == first_end

    for lv in hdr["levels"]:
        for c in lv["chunks"]:
            for name, spec in c["columns"].items():
                inside = base + spec["offset"] < first_end
                assert inside == (lv["name"] == "noyau" and name in ("copy", "position"))


def test_the_preamble_says_everything_a_prefix_needs(tmp_path):
    rep = write(tmp_path / "a.g3d", [chain_level()])
    blob = (tmp_path / "a.g3d").read_bytes()
    assert blob[:8] == MAGIC
    version, hlen, first_end = struct.unpack_from("<III", blob, 8)
    assert version == 1
    assert hlen == rep["header"]
    assert first_end == rep["first_frame_end"]
    hdr, _, _ = read_header(blob[:first_end])       # un préfixe suffit
    assert hdr["first"] == "noyau"


def test_a_prefix_shorter_than_the_header_is_refused_not_misread(tmp_path):
    write(tmp_path / "a.g3d", [chain_level()])
    blob = (tmp_path / "a.g3d").read_bytes()
    with pytest.raises(ValueError, match="trop court"):
        read_header(blob[:PREAMBLE + 10])


# --------------------------------------------------------------------------
# Intégrité
# --------------------------------------------------------------------------


def test_one_flipped_byte_in_a_column_is_detected(tmp_path):
    rep = write(tmp_path / "a.g3d", [chain_level()])
    blob = bytearray((tmp_path / "a.g3d").read_bytes())
    blob[rep["first_frame_end"] - 3] ^= 0x01
    (tmp_path / "b.g3d").write_bytes(bytes(blob))
    with pytest.raises(CorruptError):
        read(tmp_path / "b.g3d")


def test_one_flipped_byte_in_the_header_is_detected(tmp_path):
    write(tmp_path / "a.g3d", [chain_level()])
    blob = bytearray((tmp_path / "a.g3d").read_bytes())
    blob[PREAMBLE + 5] ^= 0x01
    with pytest.raises(CorruptError, match="en-tête"):
        read_header(bytes(blob))


def test_a_file_that_is_not_g3d_is_refused():
    with pytest.raises(ValueError, match="signature"):
        read_header(b"\x89PNG\r\n\x1a\n" + bytes(80))


# --------------------------------------------------------------------------
# Octree et index
# --------------------------------------------------------------------------


def test_octree_leaves_partition_the_points_and_respect_the_cap():
    rng = np.random.default_rng(2)
    pts = rng.normal(0, 1, (20_000, 3))
    leaves = octree(pts, 1_000)
    allidx = np.sort(np.concatenate(leaves))
    assert np.array_equal(allidx, np.arange(len(pts)))
    assert max(len(lf) for lf in leaves) <= 1_000


def test_a_dense_core_gives_more_leaves_than_a_sparse_shell():
    """Le découpage suit la densité, pas une grille : c'est ce qui évite des chunks vides."""
    rng = np.random.default_rng(3)
    core = rng.normal(0, 0.1, (9_000, 3))
    shell = rng.normal(0, 1, (1_000, 3))
    leaves = octree(np.concatenate([core, shell]), 500)
    near = sum(1 for lf in leaves if np.linalg.norm(np.concatenate([core, shell])[lf], axis=1).mean() < 0.5)
    assert near > len(leaves) / 2


def test_several_chunks_round_trip(tmp_path):
    lv = chain_level(n_per_copy=(2000, 1500, 1800))
    rep = write(tmp_path / "a.g3d", [lv], max_beads=700)
    assert rep["chunks"]["noyau"] > 3
    got = read(rep["path"]).levels["noyau"]
    o1, o2 = by_locus(lv.copy, lv.start), by_locus(got["copy"], got["start"])
    assert np.array_equal(got["start"][o2], lv.start[o1])


def test_chunk_spheres_contain_their_beads(tmp_path):
    lv = chain_level(n_per_copy=(2000, 1500))
    rep = write(tmp_path / "a.g3d", [lv], max_beads=600)
    got = read(rep["path"]).levels["noyau"]
    for c in got["meta"]["chunks"]:
        p = got["positions_nm"][c["first"]:c["first"] + c["n"]]
        cx, cy, cz, r = c["sphere"]
        d = np.linalg.norm(p - np.array([cx, cy, cz]), axis=1)
        assert d.max() <= r + 2 * got["meta"]["quant"]["max_error_nm"] + 0.02


def test_the_genomic_index_finds_every_bead_by_locus(tmp_path):
    lv = chain_level(n_per_copy=(500, 400), seed=4)
    rep = write(tmp_path / "a.g3d", [lv], max_beads=200)
    got = read(rep["path"]).levels["noyau"]
    idx = got["index"]
    offs = got["meta"]["index"]["copy_offsets"]
    for c in range(len(lv.copies)):
        ids = idx["ids"][offs[c]:offs[c + 1]]
        assert (got["copy"][ids] == c).all()
        starts = idx["starts"][offs[c]:offs[c + 1]]
        assert (np.diff(starts.astype(np.int64)) > 0).all()     # trié : dichotomie possible
        assert np.array_equal(got["start"][ids], starts)


# --------------------------------------------------------------------------
# Validation des entrées
# --------------------------------------------------------------------------


def test_a_bead_ending_before_it_starts_is_refused(tmp_path):
    lv = chain_level()
    lv.end = lv.end.copy()
    lv.end[7] = lv.start[7]
    with pytest.raises(ValueError, match="fin avant"):
        write(tmp_path / "a.g3d", [lv])


def test_duplicate_level_names_are_refused(tmp_path):
    with pytest.raises(ValueError, match="double"):
        write(tmp_path / "a.g3d", [chain_level(), chain_level()])


def test_the_header_carries_provenance_verbatim(tmp_path):
    hdr = {"assembly": "GRCh38", "warnings": ["simulé"]}
    write(tmp_path / "a.g3d", [chain_level()], header=hdr)
    got = read(tmp_path / "a.g3d").header
    assert got["assembly"] == "GRCh38" and got["warnings"] == ["simulé"]
    assert json.dumps(got)                                    # sérialisable tel quel


# --------------------------------------------------------------------------
# Construction : aperçu et pose
# --------------------------------------------------------------------------


class _Beads:
    def __init__(self, copy_id, start, end, nuclear_radius=5.0):
        self.copy_id, self.start, self.end = copy_id, start, end
        self.nuclear_radius = nuclear_radius
        self.labels = [f"c{i}" for i in range(int(copy_id.max()) + 1)]


def test_coarsening_never_merges_across_two_copies():
    copy = np.array([0] * 6 + [1] * 5)
    start = np.r_[np.arange(6), np.arange(5)] * 100
    b = _Beads(copy, start, start + 100)
    rng = np.random.default_rng(5)
    coords = rng.normal(0, 1, (7, 11, 3))
    lv = Level("noyau", 100, coords[0] * 1000, copy, start, start + 100,
               np.ones(11), b.labels, [100.0, 100.0])
    ap = coarsen(lv, {"beads": b, "coords_um": coords, "medoid": 0}, k=4)
    # 6 billes → groupes 4 + 2 ; 5 billes → 4 + 1
    assert ap.copy.tolist() == [0, 0, 1, 1]
    assert ap.start.tolist() == [0, 400, 0, 400]
    assert ap.end.tolist() == [400, 600, 400, 500]
    assert np.allclose(ap.positions_nm[0], coords[0, :4].mean(axis=0) * 1000)


def test_coarsened_variability_is_recomputed_not_averaged():
    """La moyenne des écarts-types filles n'est l'écart-type de rien : deux billes qui
    bougent en sens opposés ont chacune un grand écart-type, leur centroïde aucun."""
    copy = np.zeros(2, dtype=int)
    start = np.array([0, 100])
    b = _Beads(copy, start, start + 100, nuclear_radius=1.0)
    n = 50
    t = np.linspace(-1, 1, n)
    coords = np.zeros((n, 2, 3))
    coords[:, 0, 0] = 0.5 + 0.3 * t
    coords[:, 1, 0] = 0.5 - 0.3 * t
    lv = Level("noyau", 100, coords[0] * 1000, copy, start, start + 100,
               np.ones(2), b.labels, [1.0])
    ap = coarsen(lv, {"beads": b, "coords_um": coords, "medoid": 0}, k=2)
    assert ap.variability_nm[0] == pytest.approx(0.0, abs=1e-6)


def test_pose_recovers_a_known_rigid_motion_exactly():
    rng = np.random.default_rng(6)
    local = rng.normal(0, 300, (6, 3))
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    target = local @ q.T + np.array([1000.0, -200.0, 50.0])
    m, fit = pose(target, local)
    mapped = (m @ np.c_[local, np.ones(6)].T).T[:, :3]
    assert np.allclose(mapped, target, atol=1e-6)
    assert fit["rmsd_nm"] == pytest.approx(0.0, abs=0.05)
    assert fit["scale_if_allowed"] == pytest.approx(1.0, abs=1e-3)


def test_pose_reports_a_scale_it_does_not_apply():
    """Un modèle fin deux fois trop compact : la pose rigide laisse un résidu, et dit
    qu'une échelle de 2 l'annulerait — sans l'appliquer."""
    rng = np.random.default_rng(7)
    target = rng.normal(0, 400, (6, 3))
    local = target / 2.0
    m, fit = pose(target, local)
    assert fit["scale_if_allowed"] == pytest.approx(2.0, rel=1e-6)
    assert fit["rmsd_if_scaled_nm"] == pytest.approx(0.0, abs=0.05)
    assert fit["rmsd_nm"] > 50
    assert np.allclose(np.abs(np.linalg.det(m[:3, :3])), 1.0)   # rigide : aucune échelle
