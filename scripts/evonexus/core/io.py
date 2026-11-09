"""Tabular and JSON IO with checksums and a parquet/CSV fallback.

Parquet is preferred because it preserves dtypes, but a laptop install without
``pyarrow`` must still be able to run the whole pipeline. ``write_table`` picks
the best available format and always records which one it used, so the run
record states the real on-disk representation rather than an assumed one.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

try:  # pragma: no cover - depends on the install
    import pyarrow  # noqa: F401

    PARQUET_AVAILABLE = True
except Exception:  # pragma: no cover
    PARQUET_AVAILABLE = False

TABLE_SUFFIX = ".parquet" if PARQUET_AVAILABLE else ".csv.gz"
_FALLBACK_SUFFIXES = (".parquet", ".csv.gz", ".csv", ".tsv", ".tsv.gz")


def _strip_table_suffix(path: Path) -> Path:
    name = path.name
    for suffix in _FALLBACK_SUFFIXES:
        if name.endswith(suffix):
            return path.with_name(name[: -len(suffix)])
    return path


def table_path(path: str | Path) -> Path:
    """Return the concrete on-disk path for a logical table name."""
    base = _strip_table_suffix(Path(path))
    return base.with_name(base.name + TABLE_SUFFIX)


def write_table(df: pd.DataFrame, path: str | Path) -> Path:
    """Write a dataframe, creating parent directories.

    The logical ``path`` may be given with or without an extension; the
    extension actually used depends on whether parquet support is installed.
    """
    target = table_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.suffix == ".parquet":
        df.to_parquet(target, index=False)
    else:
        df.to_csv(target, index=False)
    return target


def resolve_table(path: str | Path) -> Path:
    """Find an existing table file for a logical name, trying known suffixes."""
    candidate = Path(path)
    if candidate.exists() and candidate.is_file():
        return candidate
    base = _strip_table_suffix(candidate)
    for suffix in (TABLE_SUFFIX,) + _FALLBACK_SUFFIXES:
        trial = base.with_name(base.name + suffix)
        if trial.exists():
            return trial
    raise FileNotFoundError(f"no table found for {path} (tried {_FALLBACK_SUFFIXES})")


def read_table(path: str | Path, **kwargs: Any) -> pd.DataFrame:
    """Read a table written by :func:`write_table` or supplied by a source."""
    target = resolve_table(path)
    name = target.name
    if name.endswith(".parquet"):
        return pd.read_parquet(target, **kwargs)
    if name.endswith((".tsv", ".tsv.gz")):
        return pd.read_csv(target, sep="\t", **kwargs)
    return pd.read_csv(target, **kwargs)


def read_any_table(path: str | Path, nrows: int | None = None) -> pd.DataFrame:
    """Best-effort read of an arbitrary source table (csv/tsv/parquet/xlsx)."""
    target = Path(path)
    name = target.name.lower()
    if name.endswith(".parquet"):
        df = pd.read_parquet(target)
        return df.head(nrows) if nrows else df
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(target, nrows=nrows)
    sep = "\t" if name.endswith((".tsv", ".tsv.gz", ".tab")) else None
    return pd.read_csv(target, sep=sep, engine="python", nrows=nrows)


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


def write_text(text: str, path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    return target


def sha256_file(path: str | Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_obj(obj: Any) -> str:
    """Stable checksum of a JSON-serializable object."""
    payload = json.dumps(_sanitize(obj), sort_keys=True, cls=_JSONEncoder).encode("utf-8")
    return sha256_bytes(payload)


def iter_files(root: str | Path, exclude: Iterable[str] = ()) -> list[Path]:
    excluded = tuple(exclude)
    out: list[Path] = []
    for p in sorted(Path(root).rglob("*")):
        if not p.is_file():
            continue
        if any(part in excluded for part in p.parts):
            continue
        out.append(p)
    return out
