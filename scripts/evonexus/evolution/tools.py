"""Thin wrappers around MAFFT, IQ-TREE 2 and HyPhy.

Every wrapper returns a structured result that distinguishes three outcomes:
``ok``, ``tool_missing`` and ``failed``. A missing tool is a documented
exclusion, never a silent skip and never a substituted approximation.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any


def tool_available(name: str) -> bool:
    return shutil.which(name) is not None


def _run(cmd: Sequence[str], cwd: str | Path | None = None,
         timeout: int = 7200) -> dict[str, Any]:
    try:
        proc = subprocess.run(
            list(cmd), capture_output=True, text=True, timeout=timeout,
            cwd=str(cwd) if cwd else None,
        )
    except FileNotFoundError:
        return {"status": "tool_missing", "command": list(cmd)}
    except subprocess.TimeoutExpired:
        return {"status": "failed", "command": list(cmd), "error": f"timed out after {timeout}s"}
    return {
        "status": "ok" if proc.returncode == 0 else "failed",
        "command": list(cmd),
        "returncode": proc.returncode,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def run_mafft(input_fasta: str | Path, output_fasta: str | Path,
              args: Sequence[str] = ("--auto",)) -> dict[str, Any]:
    """Align with MAFFT. MAFFT writes to stdout, so capture it and persist it."""
    if not tool_available("mafft"):
        return {"status": "tool_missing", "tool": "mafft",
                "remedy": "conda install -c bioconda mafft"}
    out = Path(output_fasta)
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        proc = subprocess.run(["mafft", *args, str(input_fasta)],
                              capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        return {"status": "failed", "tool": "mafft", "error": "timed out"}
    if proc.returncode != 0:
        return {"status": "failed", "tool": "mafft", "stderr_tail": proc.stderr[-4000:]}
    out.write_text(proc.stdout, encoding="utf-8")
    return {"status": "ok", "tool": "mafft", "output": str(out),
            "stderr_tail": proc.stderr[-2000:]}


def run_iqtree(alignment: str | Path, prefix: str | Path, model: str = "MFP",
               bootstrap: int = 1000, extra: Sequence[str] = ()) -> dict[str, Any]:
    binary = "iqtree2" if tool_available("iqtree2") else ("iqtree" if tool_available("iqtree") else None)
    if binary is None:
        return {"status": "tool_missing", "tool": "iqtree2",
                "remedy": "conda install -c bioconda iqtree"}
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "-s", str(alignment), "-m", model, "--prefix", str(prefix), *extra]
    if bootstrap:
        cmd += ["-B", str(bootstrap)]
    result = _run(cmd)
    treefile = prefix.with_suffix(".treefile")
    result["tool"] = binary
    result["treefile"] = str(treefile) if treefile.exists() else None
    return result


def run_hyphy(method: str, alignment: str | Path, tree: str | Path,
              output_json: str | Path, extra: Sequence[str] = ()) -> dict[str, Any]:
    if not tool_available("hyphy"):
        return {"status": "tool_missing", "tool": "hyphy", "method": method,
                "remedy": "conda install -c bioconda hyphy"}
    out = Path(output_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "hyphy", method.lower(),
        "--alignment", str(alignment),
        "--tree", str(tree),
        "--output", str(out),
        *extra,
    ]
    result = _run(cmd)
    result["tool"] = "hyphy"
    result["method"] = method
    result["output"] = str(out) if out.exists() else None
    return result


def parse_hyphy_sites(json_path: str | Path, method: str,
                      p_threshold: float = 0.05) -> dict[str, Any]:
    """Extract a site-level summary from a HyPhy result JSON."""
    path = Path(json_path)
    if not path.exists():
        return {"status": "missing", "method": method}
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"status": "unparsable", "method": method, "error": str(exc)}

    content = (doc.get("MLE") or {}).get("content") or {}
    headers = [h[0] for h in (doc.get("MLE") or {}).get("headers", [])]
    rows = content.get("0") or []
    if not rows:
        return {"status": "no_sites", "method": method,
                "test_result": doc.get("test results")}

    p_index = next((i for i, h in enumerate(headers) if "p-value" in h.lower()), None)
    n_sig = 0
    if p_index is not None:
        n_sig = sum(1 for row in rows if len(row) > p_index and row[p_index] <= p_threshold)
    return {
        "status": "ok",
        "method": method,
        "n_sites": len(rows),
        "n_significant": n_sig,
        "p_threshold": p_threshold,
        "headers": headers,
        "test_result": doc.get("test results"),
    }
