"""Inventory of whatever a source deposit actually contains.

The inventory is deliberately opinion-free: it reports real headers, real row
counts and a guess at each table's role, and it never renames or reshapes
anything. A human reads the inventory, fills in ``config/metadata_mapping.yaml``
and marks it approved. Until then harmonization refuses to run.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from evonexus.core.io import read_any_table, sha256_file, write_json

TABULAR_SUFFIXES = {".csv", ".tsv", ".txt", ".tab", ".xlsx", ".xls", ".parquet"}
SEQUENCE_SUFFIXES = {".fa", ".fasta", ".fna", ".faa", ".fas", ".gb", ".gbk", ".genbank"}
ALIGNMENT_SUFFIXES = {".aln", ".phy", ".nex", ".nexus", ".msa", ".afa"}
TREE_SUFFIXES = {".nwk", ".newick", ".tree", ".treefile", ".nex"}
ARCHIVE_SUFFIXES = {".zip", ".gz", ".tar", ".tgz", ".bz2", ".xz"}

_ROLE_HINTS = {
    "detections": [r"detect", r"abundance", r"rpkm", r"count.*matrix", r"votu.*table", r"presence"],
    "libraries": [r"librar", r"sample.*meta", r"runinfo", r"sequenc.*meta", r"metadata.*sample"],
    "hosts": [r"host", r"\bbat\b", r"specimen", r"pool", r"animal"],
    "sites": [r"site", r"location", r"geograph", r"sampling.*site"],
    "viruses": [r"virus", r"viral", r"votu", r"taxonom", r"catalog", r"lineage"],
    "ecology": [r"ecolog", r"environment", r"climate", r"habitat", r"land.*cover"],
}


def _classify_kind(path: Path) -> str:
    suffixes = [s.lower() for s in path.suffixes]
    last = suffixes[-1] if suffixes else ""
    stem_suffix = suffixes[-2] if len(suffixes) > 1 else ""
    effective = stem_suffix if last in ARCHIVE_SUFFIXES and stem_suffix else last
    if effective in TABULAR_SUFFIXES:
        return "table"
    if effective in SEQUENCE_SUFFIXES:
        return "sequences"
    if effective in ALIGNMENT_SUFFIXES:
        return "alignment"
    if effective in TREE_SUFFIXES:
        return "tree"
    if last in ARCHIVE_SUFFIXES:
        return "archive"
    return "other"


def _guess_role(path: Path, columns: List[str]) -> Dict[str, Any]:
    haystack = (path.name + " " + " ".join(columns)).lower()
    scores: Dict[str, int] = {}
    for role, patterns in _ROLE_HINTS.items():
        scores[role] = sum(1 for pat in patterns if re.search(pat, haystack))
    best = max(scores, key=lambda k: scores[k])
    return {
        "guess": best if scores[best] > 0 else "unclassified",
        "confidence": "low" if scores[best] <= 1 else "medium",
        "scores": scores,
        "note": "A guess only. Confirm it in config/metadata_mapping.yaml before use.",
    }


def describe_table(path: Path, preview_rows: int = 5) -> Dict[str, Any]:
    entry: Dict[str, Any] = {"readable": False}
    try:
        head = read_any_table(path, nrows=2000)
    except Exception as exc:
        entry["error"] = f"{type(exc).__name__}: {exc}"
        return entry
    entry["readable"] = True
    entry["columns"] = [str(c) for c in head.columns]
    entry["n_columns"] = int(head.shape[1])
    entry["rows_sampled"] = int(head.shape[0])
    entry["dtypes"] = {str(c): str(t) for c, t in head.dtypes.items()}
    entry["null_fraction"] = {
        str(c): round(float(head[c].isna().mean()), 4) for c in head.columns
    }
    entry["example_values"] = {
        str(c): [None if pd.isna(v) else str(v) for v in head[c].dropna().head(3).tolist()]
        for c in head.columns
    }
    entry["preview"] = head.head(preview_rows).astype(str).to_dict(orient="records")
    entry["role_hint"] = _guess_role(path, entry["columns"])
    return entry


def inventory_directory(root: str | Path, checksum: bool = True,
                        max_checksum_bytes: int = 5 * (1 << 30)) -> Dict[str, Any]:
    """Walk a raw-data directory and describe everything found."""
    base = Path(root)
    files: List[Dict[str, Any]] = []
    for path in sorted(base.rglob("*")):
        if not path.is_file() or path.name.startswith("."):
            continue
        kind = _classify_kind(path)
        size = path.stat().st_size
        entry: Dict[str, Any] = {
            "path": str(path.relative_to(base)),
            "absolute_path": str(path),
            "kind": kind,
            "bytes": size,
        }
        if checksum and size <= max_checksum_bytes:
            entry["sha256"] = sha256_file(path)
        if kind == "table":
            entry.update(describe_table(path))
        files.append(entry)

    by_kind: Dict[str, int] = {}
    for f in files:
        by_kind[f["kind"]] = by_kind.get(f["kind"], 0) + 1

    tables = [f for f in files if f["kind"] == "table"]
    roles: Dict[str, List[str]] = {}
    for t in tables:
        role = (t.get("role_hint") or {}).get("guess", "unclassified")
        roles.setdefault(role, []).append(t["path"])

    return {
        "root": str(base),
        "n_files": len(files),
        "total_bytes": sum(f["bytes"] for f in files),
        "by_kind": by_kind,
        "candidate_roles": roles,
        "files": files,
        "expectations": _expectations(roles, by_kind),
    }


def _expectations(roles: Dict[str, List[str]], by_kind: Dict[str, int]) -> Dict[str, Any]:
    """State plainly which pipeline prerequisites the deposit appears to meet."""
    return {
        "library_level_metadata": bool(roles.get("libraries")),
        "host_metadata": bool(roles.get("hosts")),
        "site_metadata": bool(roles.get("sites")),
        "viral_catalog": bool(roles.get("viruses")),
        "library_by_lineage_detections": bool(roles.get("detections")),
        "alignments": by_kind.get("alignment", 0) > 0,
        "trees": by_kind.get("tree", 0) > 0,
        "verdict": (
            "processed_matrix_candidate"
            if roles.get("detections") and roles.get("libraries")
            else "raw_read_reconstruction_may_be_required"
        ),
        "caveat": (
            "Presence of a detection-like table does not prove that every "
            "library was screened against a common reference panel. Gate G3 "
            "checks that separately."
        ),
    }


def write_inventory(inventory: Dict[str, Any], out_dir: str | Path) -> Dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    json_path = write_json(inventory, out / "inventory.json")

    lines = [
        "# Source inventory",
        "",
        f"Root: `{inventory['root']}`",
        f"Files: {inventory['n_files']}  ·  Total size: {inventory['total_bytes'] / 1e6:.1f} MB",
        "",
        "## File kinds",
        "",
        "| Kind | Count |",
        "| --- | --- |",
    ]
    for kind, n in sorted(inventory["by_kind"].items()):
        lines.append(f"| {kind} | {n} |")

    lines += ["", "## Candidate roles (guesses — confirm before use)", ""]
    for role, paths in sorted(inventory["candidate_roles"].items()):
        lines.append(f"- **{role}**: " + ", ".join(f"`{p}`" for p in paths))

    lines += ["", "## Table headers", ""]
    for f in inventory["files"]:
        if f["kind"] != "table":
            continue
        lines.append(f"### `{f['path']}`")
        lines.append("")
        if not f.get("readable"):
            lines.append(f"Not readable: {f.get('error', 'unknown error')}")
            lines.append("")
            continue
        lines.append(f"{f['n_columns']} columns, {f['rows_sampled']} rows sampled.")
        lines.append("")
        lines.append("| Column | dtype | null fraction | example |")
        lines.append("| --- | --- | --- | --- |")
        for col in f["columns"]:
            ex = ", ".join(f["example_values"].get(col, [])[:2]) or "—"
            lines.append(
                f"| `{col}` | {f['dtypes'].get(col, '?')} | "
                f"{f['null_fraction'].get(col, 0):.2f} | {ex} |"
            )
        lines.append("")

    lines += ["", "## Prerequisite expectations", ""]
    for key, value in inventory["expectations"].items():
        lines.append(f"- **{key}**: {value}")
    lines.append("")

    md_path = out / "inventory.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    return {"json": json_path, "markdown": md_path}
