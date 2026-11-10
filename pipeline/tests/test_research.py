"""Mode Recherche (semaine 17) : la session sur un petit `.mcool` planté, et le service local."""

from __future__ import annotations

import json
import threading
import urllib.error
import urllib.request

import numpy as np
import pytest

pytest.importorskip("cooler")

from geno_pipeline.research import planted
from geno_pipeline.research import session as S
from geno_pipeline.research.serve import serve


@pytest.fixture(scope="module")
def mcool(tmp_path_factory):
    import logging
    import warnings

    logging.disable(logging.INFO)
    warnings.simplefilter("ignore")
    out = tmp_path_factory.mktemp("plante") / "p.mcool"
    planted.plant(out, n=(160, 120), total=6_000_000, factors=(1, 2))
    return out


def test_inspect_lists_what_cooler_wrote(mcool):
    info = S.inspect(mcool)
    assert info["chroms"] == [["chrP1", 160 * 50_000], ["chrP2", 120 * 50_000]]
    assert [(r["resolution"], r["bins"], r["balanced"]) for r in info["resolutions"]] == [(50_000, 280, True), (100_000, 140, True)]


def test_the_structure_is_the_planted_one_and_not_a_shuffle(mcool):
    r = S.reconstruct(mcool, S.Params(50_000, "chrP1", balance="aucun"))
    t = planted.truth(mcool, "chrP1")[r.start // 50_000]
    assert S.compare(r.positions_nm, t) < 0.25
    assert S.compare(r.positions_nm, t[np.random.default_rng(0).permutation(len(t))]) > 0.8
    # Convention d'échelle : deux voisines se touchent (médiane).
    step = np.linalg.norm(np.diff(r.positions_nm, axis=0), axis=1)
    assert abs(np.median(step) / (2 * r.stats["radius_nm"]) - 1) < 0.05


def test_orientation_is_a_convention_not_lapack_s_mood():
    x = np.random.default_rng(1).normal(size=(50, 3))
    a = S.orient(x.copy())
    b = S.orient(x.copy() * np.array([-1, 1, -1]))
    assert np.allclose(a, b) and np.allclose(a.mean(axis=0), 0)


def test_manifest_reproduces_bytes_and_refuses_another_file(mcool, tmp_path):
    m = S.run(mcool, S.Params(100_000, "chrP2", 0.4, "ice"), tmp_path / "a")
    r = S.reproduce(json.loads(json.dumps(m)), mcool, tmp_path / "b")
    assert r["identical"] and r["versions_differ"] == {}
    assert (tmp_path / "a" / "structure.g3d").read_bytes() == (tmp_path / "b" / "structure.g3d").read_bytes()
    # Témoin : un paramètre qui bouge change les octets.
    m2 = json.loads(json.dumps(m))
    m2["params"]["alpha"] = 0.41
    assert not S.reproduce(m2, mcool, tmp_path / "c")["identical"]
    other = tmp_path / "other.mcool"
    planted.plant(other, n=(40, 40), total=400_000, factors=(1,), seed=5)
    with pytest.raises(S.SessionError, match="pas celui du manifeste"):
        S.reproduce(m, other, tmp_path / "d")
    with pytest.raises(S.SessionError):
        S.reproduce({"format": "autre"}, mcool, tmp_path / "e")


def test_refusals_say_why(mcool, tmp_path):
    for p, word in ((S.Params(25_000, "chrP1"), "absente"), (S.Params(50_000, "chr9"), "inconnue"),
                    (S.Params(50_000, "chrP1", balance="magie"), "équilibrage"),
                    (S.Params(50_000, "chrP1:1-200,000"), "trop peu")):
        with pytest.raises(S.SessionError, match=word):
            S.reconstruct(mcool, p)
    bad = tmp_path / "x.mcool"
    bad.write_bytes(b"pas du hdf5")
    with pytest.raises(S.SessionError):
        S.inspect(bad)


def test_service_runs_a_session_and_refuses_foreign_origins(mcool, tmp_path):
    srv = serve(tmp_path / "work", None, 0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"

    def call(path, data=None, headers=None, raw=False):
        req = urllib.request.Request(base + path, data=data, headers=headers or {}, method="POST" if data is not None else "GET")
        with urllib.request.urlopen(req) as r:
            b = r.read()
            return b if raw else json.loads(b)

    try:
        up = call("/api/upload?name=p.mcool", mcool.read_bytes())
        assert up["inspect"]["resolutions"][0]["resolution"] == 50_000
        job = call("/api/run", json.dumps({"file": up["file"], "params": {"resolution": 100_000, "region": "chrP1", "balance": "aucun"}}).encode())["job"]
        for _ in range(600):
            st = call(f"/api/job/{job}")
            if st["state"] not in ("en attente", "en cours"):
                break
            threading.Event().wait(0.05)
        assert st["state"] == "fait", st
        g3d = call(f"/sessions/{job}/structure.g3d", raw=True)
        import hashlib

        assert hashlib.sha256(g3d).hexdigest() == st["manifest"]["output"]["sha256"]
        rng = urllib.request.Request(base + f"/sessions/{job}/structure.g3d", headers={"Range": "bytes=0-63"})
        with urllib.request.urlopen(rng) as r:
            assert r.status == 206 and r.read() == g3d[:64]
        with pytest.raises(urllib.error.HTTPError) as e:
            call("/api/open", json.dumps({"path": str(mcool)}).encode(), {"Origin": "https://ailleurs.example"})
        assert e.value.code == 403
        with pytest.raises(urllib.error.HTTPError) as e:
            call("/api/run", json.dumps({"file": "inconnu", "params": {}}).encode())
        assert e.value.code == 400
    finally:
        srv.shutdown()
