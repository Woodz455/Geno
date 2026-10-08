import json
from pathlib import Path

import pytest

from geno_pipeline.intervals import Store, write_store
from geno_pipeline.parsers import read_track

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture(scope="session")
def store(tmp_path_factory) -> Store:
    spec = json.loads((FIXTURES / "tracks.json").read_text(encoding="utf-8"))
    out = tmp_path_factory.mktemp("store")
    write_store(
        out,
        {t["name"]: read_track(t, FIXTURES) for t in spec["tracks"]},
        {"assembly": spec["assembly"], "source": spec["source"]},
    )
    s = Store(out)
    yield s
    s.close()
