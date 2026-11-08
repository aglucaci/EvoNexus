"""Configuration loading with explicit, auditable merge semantics."""

from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Dict, Mapping

import yaml

CONFIG_DIR_ENV = "EVONEXUS_CONFIG_DIR"
_DEFAULT_CONFIG_DIR = Path(__file__).resolve().parents[3] / "config"


def config_dir() -> Path:
    """Directory holding the YAML configuration files."""
    return Path(os.environ.get(CONFIG_DIR_ENV, _DEFAULT_CONFIG_DIR))


def load_yaml(path: str | Path) -> Dict[str, Any]:
    """Load a YAML document, returning ``{}`` for an empty file."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"configuration file not found: {p}")
    with p.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return data or {}


def load_config(name: str) -> Dict[str, Any]:
    """Load a named configuration file from the configuration directory.

    ``name`` may be given with or without the ``.yaml`` suffix, or as an
    explicit path.
    """
    candidate = Path(name)
    if candidate.suffix in {".yaml", ".yml"} and candidate.exists():
        return load_yaml(candidate)
    stem = candidate.stem if candidate.suffix in {".yaml", ".yml"} else name
    return load_yaml(config_dir() / f"{stem}.yaml")


def deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``override`` into ``base`` without mutating either."""
    out: Dict[str, Any] = copy.deepcopy(dict(base))
    for key, value in override.items():
        if key in out and isinstance(out[key], Mapping) and isinstance(value, Mapping):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def resolve_variant(config: Mapping[str, Any], variant: str) -> Dict[str, Any]:
    """Resolve a named variant from a ``variants:`` block.

    A variant may declare ``inherits: <other>`` plus an ``overrides:`` mapping.
    Inheritance chains are followed and cycles are rejected.
    """
    variants = config.get("variants") or {}
    if variant not in variants:
        available = ", ".join(sorted(variants)) or "<none>"
        raise KeyError(f"unknown variant {variant!r}; available: {available}")

    chain: list[str] = []
    cursor = variant
    while cursor is not None:
        if cursor in chain:
            raise ValueError(f"cyclic variant inheritance: {' -> '.join(chain + [cursor])}")
        chain.append(cursor)
        cursor = (variants.get(cursor) or {}).get("inherits")

    resolved: Dict[str, Any] = {}
    for name in reversed(chain):
        spec = dict(variants.get(name) or {})
        spec.pop("inherits", None)
        overrides = spec.pop("overrides", None) or {}
        spec.pop("description", None)
        resolved = deep_merge(resolved, spec)
        resolved = deep_merge(resolved, overrides)
    resolved["_variant"] = variant
    resolved["_inheritance_chain"] = list(reversed(chain))
    return resolved
