"""Resolve source table columns onto the canonical schema.

The resolver never guesses silently. It produces a report listing, for every
canonical field, which source column was used and why — or that nothing was
found. A mapping that leaves a required field unresolved is an error, not a
column full of nulls.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from evonexus.core.io import read_json, write_text

CANONICAL_ENTITIES = ("sites", "hosts", "libraries", "viruses", "detections", "ecology")


class MappingError(RuntimeError):
    """Raised when the source tables cannot be mapped onto the schema."""


def _normalize(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


@dataclass
class FieldResolution:
    canonical: str
    source_column: Optional[str]
    required: bool
    method: str              # exact | normalized | absent
    candidates: List[str] = field(default_factory=list)


@dataclass
class EntityResolution:
    entity: str
    source_table: Optional[str]
    source_path: Optional[str]
    fields: List[FieldResolution] = field(default_factory=list)
    table_selection_method: str = "unresolved"

    @property
    def missing_required(self) -> List[str]:
        return [f.canonical for f in self.fields if f.required and f.source_column is None]

    def column_map(self) -> Dict[str, str]:
        return {f.source_column: f.canonical for f in self.fields if f.source_column}


@dataclass
class ResolvedMapping:
    entities: Dict[str, EntityResolution]
    approved: bool
    approved_by: Optional[str]
    vocabularies: Dict[str, Any]
    missing_tokens: List[str]
    extras: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    @property
    def blocking_problems(self) -> List[str]:
        problems: List[str] = []
        for name, ent in self.entities.items():
            if ent.source_path is None:
                problems.append(f"{name}: no source table identified")
                continue
            for miss in ent.missing_required:
                problems.append(f"{name}: required field '{miss}' has no source column")
        return problems


def _table_entries(inventory: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [f for f in inventory.get("files", [])
            if f.get("kind") == "table" and f.get("readable")]


def _select_table(entity: str, declared: Optional[str], tables: List[Dict[str, Any]],
                  overrides: Dict[str, Any], allow_role_guess: bool) -> tuple[Optional[Dict[str, Any]], str]:
    override = ((overrides.get("entities") or {}).get(entity) or {}).get("source_table")
    wanted = declared or override
    if wanted:
        for t in tables:
            if Path(t["path"]).name == Path(wanted).name or t["path"] == wanted:
                return t, "declared" if declared else "override"
        return None, f"declared table {wanted!r} not present in the inventory"
    if allow_role_guess:
        scored = [
            t for t in tables
            if (t.get("role_hint") or {}).get("guess") == entity
        ]
        if len(scored) == 1:
            return scored[0], "role_guess"
        if len(scored) > 1:
            return None, f"{len(scored)} tables guessed as {entity!r}; declare one explicitly"
    return None, "no source table declared"


def _resolve_fields(spec: Dict[str, Any], columns: List[str]) -> List[FieldResolution]:
    by_exact = {c: c for c in columns}
    by_norm: Dict[str, str] = {}
    for c in columns:
        by_norm.setdefault(_normalize(c), c)

    out: List[FieldResolution] = []
    for canonical, rule in (spec.get("fields") or {}).items():
        candidates = list(rule.get("candidates") or [])
        required = bool(rule.get("required", False))
        chosen, method = None, "absent"
        for cand in candidates:
            if cand in by_exact:
                chosen, method = by_exact[cand], "exact"
                break
        if chosen is None:
            for cand in candidates:
                key = _normalize(cand)
                if key in by_norm:
                    chosen, method = by_norm[key], "normalized"
                    break
        out.append(FieldResolution(canonical, chosen, required, method, candidates))
    return out


def resolve_mapping(
    inventory: Dict[str, Any],
    mapping_config: Dict[str, Any],
    overrides_path: str | Path | None = None,
    allow_role_guess: bool = False,
    require_approval: bool = True,
) -> ResolvedMapping:
    """Resolve the configured mapping against a concrete inventory."""
    overrides: Dict[str, Any] = {}
    if overrides_path and Path(overrides_path).exists():
        overrides = read_json(overrides_path)

    approved = bool(mapping_config.get("approved") or overrides.get("approved"))
    approved_by = mapping_config.get("approved_by") or overrides.get("approved_by")

    tables = _table_entries(inventory)
    entities: Dict[str, EntityResolution] = {}
    warnings: List[str] = []

    for entity in CANONICAL_ENTITIES:
        spec = (mapping_config.get("entities") or {}).get(entity) or {}
        table, method = _select_table(
            entity, spec.get("source_table"), tables, overrides, allow_role_guess
        )
        res = EntityResolution(entity=entity, source_table=None, source_path=None,
                               table_selection_method=method)
        if table is not None:
            res.source_table = table["path"]
            res.source_path = table.get("absolute_path") or table["path"]
            res.fields = _resolve_fields(spec, table.get("columns", []))
        entities[entity] = res

    if require_approval and not approved:
        warnings.append(
            "config/metadata_mapping.yaml is not approved. Review "
            "results/qc/inventory.md, complete the mapping and set "
            "`approved: true` before harmonizing real data."
        )

    return ResolvedMapping(
        entities=entities,
        approved=approved,
        approved_by=approved_by,
        vocabularies=mapping_config.get("vocabularies") or {},
        missing_tokens=list(mapping_config.get("missing_tokens") or []),
        extras=overrides.get("extras") or {},
        warnings=warnings,
    )


def write_mapping_report(mapping: ResolvedMapping, path: str | Path) -> Path:
    lines = [
        "# Metadata mapping report",
        "",
        f"Mapping approved: **{mapping.approved}**"
        + (f" (by {mapping.approved_by})" if mapping.approved_by else ""),
        "",
    ]
    if mapping.warnings:
        lines += ["## Warnings", ""] + [f"- {w}" for w in mapping.warnings] + [""]

    problems = mapping.blocking_problems
    lines += ["## Blocking problems", ""]
    lines += [f"- {p}" for p in problems] if problems else ["None."]
    lines.append("")

    for name, ent in mapping.entities.items():
        lines += [f"## `{name}`", ""]
        lines.append(f"- Source table: `{ent.source_table or 'NOT RESOLVED'}`")
        lines.append(f"- Selection method: {ent.table_selection_method}")
        lines.append("")
        if not ent.fields:
            lines += ["No fields resolved.", ""]
            continue
        lines += [
            "| Canonical field | Source column | Required | Match |",
            "| --- | --- | --- | --- |",
        ]
        for f in ent.fields:
            lines.append(
                f"| `{f.canonical}` | "
                f"{'`' + f.source_column + '`' if f.source_column else '**—**'} | "
                f"{'yes' if f.required else 'no'} | {f.method} |"
            )
        lines.append("")

    return write_text("\n".join(lines), path)
