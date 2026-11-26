"""Standalone checks for the installed package.

This suite exists to be run against an installed distribution with the
EvoNexus repository *not* on ``sys.path``. Its value is the isolation, not
the depth of its assertions: a stray ``from evonexus...`` import left behind
in the package would pass every other CI job silently and fail here.
"""

from __future__ import annotations

import json

import numpy as np

import evonexus_bayes
from evonexus_bayes._json import read_json, write_json


def test_version_is_exposed() -> None:
    assert evonexus_bayes.__version__


def test_package_does_not_import_evonexus() -> None:
    import sys

    assert not [name for name in sys.modules if name.split(".")[0] == "evonexus"]


def test_json_roundtrip_handles_numpy_and_non_finite(tmp_path) -> None:
    payload = {
        "draws": [np.int64(0), np.int64(1)],
        "scale": np.float64(1.5),
        "diverged": float("nan"),
        "bound": float("inf"),
        # A bare ndarray has `.item()`, but it raises for size > 1, so the
        # encoder falls through to `str`. Pinned here because it is the
        # pipeline's long-standing behaviour, not an intended design.
        "array": np.arange(3),
    }
    target = write_json(payload, tmp_path / "posterior.json")
    assert json.loads(target.read_text()) == {
        "array": "[0 1 2]",
        "bound": None,
        "diverged": None,
        "draws": [0, 1],
        "scale": 1.5,
    }
    assert read_json(target)["scale"] == 1.5
