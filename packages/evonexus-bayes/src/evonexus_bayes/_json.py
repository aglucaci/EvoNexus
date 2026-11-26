"""JSON writing with the same conventions EvoNexus uses.

These three functions are deliberately a copy of ``evonexus.core.io``'s JSON
half rather than an import of it: the package must stand alone, and the
alternative — depending on the pipeline for ``read_json`` and ``write_json`` —
would invert the dependency this package exists to remove.

The copy is small and the drift risk is real, so ``tests/test_models.py`` in
the EvoNexus repository asserts that both implementations produce
byte-identical output for a numpy-laden object. If you change the encoding
here, that test fails until the pipeline side changes too.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


class _JSONEncoder(json.JSONEncoder):
    def default(self, o: Any) -> Any:  # pragma: no cover - defensive
        if isinstance(o, Path):
            return str(o)
        if hasattr(o, "item"):
            try:
                return o.item()
            except Exception:
                pass
        if isinstance(o, (set, frozenset)):
            return sorted(o)
        return str(o)


def _sanitize(obj: Any) -> Any:
    """Replace non-finite floats so the JSON stays strictly valid."""
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def write_json(obj: Any, path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        json.dump(_sanitize(obj), handle, indent=2, sort_keys=True, cls=_JSONEncoder)
        handle.write("\n")
    return target


def read_json(path: str | Path) -> Any:
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)
