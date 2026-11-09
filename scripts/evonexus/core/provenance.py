"""Run records: what ran, on what inputs, with which software and seed."""

from __future__ import annotations

import platform
import subprocess
import sys
import time
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from evonexus.core.io import sha256_file, write_json

_TRACKED_PACKAGES = (
    "numpy",
    "pandas",
    "scipy",
    "pymc",
    "arviz",
    "pytensor",
    "matplotlib",
    "yaml",
    "pyarrow",
    "Bio",
    "networkx",
)


def utc_now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _git(*args: str) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args],
            capture_output=True,
            text=True,
            timeout=15,
            cwd=Path(__file__).resolve().parents[3],
        )
    except Exception:
        return None
    if out.returncode != 0:
        return None
    return out.stdout.strip() or None


def package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in _TRACKED_PACKAGES:
        try:
            module = __import__(name)
        except Exception:
            versions[name] = "not installed"
            continue
        versions[name] = str(getattr(module, "__version__", "unknown"))
    return versions


def record_environment() -> dict[str, Any]:
    """Everything needed to say which software produced an output."""
    return {
        "timestamp_utc": utc_now(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": package_versions(),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
    }


@dataclass
class RunRecord:
    """Accumulates provenance for one pipeline step and writes it as JSON."""

    step: str
    params: dict[str, Any] = field(default_factory=dict)
    inputs: list[dict[str, Any]] = field(default_factory=list)
    outputs: list[dict[str, Any]] = field(default_factory=list)
    messages: list[dict[str, str]] = field(default_factory=list)
    status: str = "running"
    synthetic: bool = False
    started_at: str = field(default_factory=utc_now)
    _t0: float = field(default_factory=time.monotonic, repr=False)

    # -- registration ----------------------------------------------------
    def add_input(self, path: str | Path, role: str = "input", checksum: bool = True) -> None:
        self.inputs.append(self._describe(path, role, checksum))

    def add_output(self, path: str | Path, role: str = "output", checksum: bool = True) -> None:
        self.outputs.append(self._describe(path, role, checksum))

    def add_inputs(self, paths: Iterable[str | Path], role: str = "input") -> None:
        for p in paths:
            self.add_input(p, role)

    @staticmethod
    def _describe(path: str | Path, role: str, checksum: bool) -> dict[str, Any]:
        p = Path(path)
        entry: dict[str, Any] = {"path": str(p), "role": role, "exists": p.exists()}
        if p.exists() and p.is_file():
            entry["bytes"] = p.stat().st_size
            if checksum and entry["bytes"] <= 2 * (1 << 30):
                entry["sha256"] = sha256_file(p)
        return entry

    # -- messages --------------------------------------------------------
    def note(self, message: str) -> None:
        self.messages.append({"level": "info", "message": message})
        print(f"[evonexus] {message}", flush=True)

    def warn(self, message: str) -> None:
        self.messages.append({"level": "warning", "message": message})
        print(f"[evonexus][WARN] {message}", file=sys.stderr, flush=True)

    def error(self, message: str) -> None:
        self.messages.append({"level": "error", "message": message})
        print(f"[evonexus][ERROR] {message}", file=sys.stderr, flush=True)

    # -- finalization ----------------------------------------------------
    def finish(self, path: str | Path, status: str = "ok") -> Path:
        self.status = status
        payload = {
            "step": self.step,
            "status": self.status,
            "synthetic": self.synthetic,
            "started_at": self.started_at,
            "finished_at": utc_now(),
            "duration_seconds": round(time.monotonic() - self._t0, 3),
            "params": self.params,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "messages": self.messages,
            "environment": record_environment(),
            "evonexus_version": __import__("evonexus").__version__,
        }
        return write_json(payload, path)
