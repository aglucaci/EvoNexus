"""Shared fixtures. A single small synthetic dataset backs most tests."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def synthetic_raw(tmp_path_factory) -> Path:
    from evonexus.simulate import simulate

    out = tmp_path_factory.mktemp("synthetic") / "raw" / "synthetic"
    simulate(out, n_sites=8, n_libraries=40, n_viruses=12, seed=7)
    return out.parent


@pytest.fixture(scope="session")
def inventory(synthetic_raw, tmp_path_factory):
    from evonexus.acquire.inventory import inventory_directory

    return inventory_directory(synthetic_raw)


@pytest.fixture(scope="session")
def canonical(synthetic_raw, inventory):
    from evonexus.core.config import load_config
    from evonexus.harmonize.canonical import build_canonical
    from evonexus.harmonize.mapping import resolve_mapping

    mapping = resolve_mapping(
        inventory,
        load_config("metadata_mapping"),
        overrides_path=synthetic_raw / "synthetic" / "mapping_overrides.json",
        allow_role_guess=True,
        require_approval=False,
    )
    return build_canonical(synthetic_raw, mapping, load_config("viral_filters"))


@pytest.fixture(scope="session")
def tables(canonical):
    return canonical["tables"]


@pytest.fixture(scope="session")
def analysis_frame(tables):
    from evonexus.core.config import load_config
    from evonexus.quantify.matrix import build_analysis_frame

    return build_analysis_frame(tables, load_config("viral_filters"))
