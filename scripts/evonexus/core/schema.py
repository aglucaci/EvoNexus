"""Canonical schemas for the harmonized datasets.

The schema layer encodes one rule that the rest of the pipeline relies on:
*unknown is not zero*. Habitat, host species and detection status each carry an
explicit ``unknown`` / ``unassessable`` state that cannot be confused with an
observed category, and validation refuses tables that lose that distinction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Sequence

import pandas as pd


class SchemaError(ValueError):
    """Raised when a table violates its canonical schema."""


@dataclass(frozen=True)
class Column:
    name: str
    dtype: str                       # string | float | int | bool | category
    required: bool = True
    nullable: bool = True
    allowed: Sequence[str] | None = None
    description: str = ""


@dataclass(frozen=True)
class TableSchema:
    name: str
    unit: str
    primary_key: Sequence[str]
    columns: Sequence[Column]
    foreign_keys: Dict[str, str] = field(default_factory=dict)

    def column(self, name: str) -> Column | None:
        for col in self.columns:
            if col.name == name:
                return col
        return None


HABITAT_LEVELS = ("settlement", "wild", "unknown")
DETECTION_LEVELS = ("positive", "negative", "unassessable")
POOL_LEVELS = ("individual", "pool", "unknown")

CANONICAL_SCHEMAS: Dict[str, TableSchema] = {
    "sites": TableSchema(
        name="sites",
        unit="sampling site",
        primary_key=("site_id",),
        columns=[
            Column("site_id", "string", nullable=False, description="Canonical site identifier"),
            Column("region", "string", required=False),
            Column("country", "string", required=False),
            Column("habitat", "category", allowed=HABITAT_LEVELS,
                   description="Settlement / wild / unknown; unknown is excluded from the contrast"),
            Column("habitat_raw", "string", required=False, description="Verbatim source value"),
            Column("latitude", "float", required=False),
            Column("longitude", "float", required=False),
            Column("geo_resolution", "string", required=False,
                   description="Stated precision of the coordinates"),
        ],
    ),
    "hosts": TableSchema(
        name="hosts",
        unit="animal or documented pool",
        primary_key=("host_id",),
        foreign_keys={"site_id": "sites.site_id"},
        columns=[
            Column("host_id", "string", nullable=False),
            Column("site_id", "string"),
            Column("host_species", "string",
                   description="Verbatim species, or 'mixed' / 'unknown'"),
            Column("host_species_resolved", "string",
                   description="Harmonized species name used for modelling"),
            Column("host_genus", "string", required=False),
            Column("host_family", "string", required=False),
            Column("pool_status", "category", allowed=POOL_LEVELS),
            Column("pool_size", "float", required=False,
                   description="Number of animals; NaN when undocumented"),
            Column("mixed_species_pool", "bool", required=False,
                   description="True when the pool spans more than one species"),
            Column("collection_date", "string", required=False),
        ],
    ),
    "libraries": TableSchema(
        name="libraries",
        unit="sequencing library",
        primary_key=("library_id",),
        foreign_keys={"host_id": "hosts.host_id", "site_id": "sites.site_id"},
        columns=[
            Column("library_id", "string", nullable=False),
            Column("host_id", "string"),
            Column("site_id", "string"),
            Column("tissue", "string"),
            Column("tissue_group", "string", description="Controlled tissue grouping"),
            Column("read_depth", "float", description="Total reads reported by the source"),
            Column("effective_depth", "float",
                   description="Depth denominator used in the model; see viral_filters.yaml"),
            Column("effective_depth_definition", "string"),
            Column("sequencing_batch", "string", required=False),
            Column("sra_run", "string", required=False),
            Column("bioproject", "string", required=False),
            Column("collection_date", "string", required=False),
            Column("collection_time_numeric", "float", required=False,
                   description="Decimal years since the first collection date"),
            Column("screened_panel", "string", required=False,
                   description="Identifier of the reference panel this library was screened against"),
            Column("is_negative_control", "bool", required=False),
        ],
    ),
    "viruses": TableSchema(
        name="viruses",
        unit="viral lineage",
        primary_key=("virus_id",),
        columns=[
            Column("virus_id", "string", nullable=False),
            Column("virus_family", "string"),
            Column("virus_name", "string", required=False),
            Column("classification_evidence", "string", required=False),
            Column("representative_accession", "string", required=False),
            Column("mammal_associated", "bool", required=False),
            Column("eligible", "bool", required=False,
                   description="Passes viral_filters.yaml lineage eligibility"),
            Column("eligibility_reason", "string", required=False),
        ],
    ),
    "detections": TableSchema(
        name="detections",
        unit="library x viral lineage",
        primary_key=("library_id", "virus_id"),
        foreign_keys={"library_id": "libraries.library_id", "virus_id": "viruses.virus_id"},
        columns=[
            Column("library_id", "string", nullable=False),
            Column("virus_id", "string", nullable=False),
            Column("detection", "category", allowed=DETECTION_LEVELS, nullable=False,
                   description="positive / negative / unassessable — never silently zero"),
            Column("y", "float", required=False,
                   description="1.0, 0.0 or NaN mirroring `detection`"),
            Column("mapped_reads", "float", required=False),
            Column("coverage_breadth", "float", required=False),
            Column("abundance", "float", required=False),
            Column("evidence", "string", required=False,
                   description="published_table | read_mapping | not_screened"),
            Column("flag", "string", required=False,
                   description="Contamination / index-hopping flags, semicolon separated"),
        ],
    ),
    "ecology": TableSchema(
        name="ecology",
        unit="site",
        primary_key=("site_id",),
        foreign_keys={"site_id": "sites.site_id"},
        columns=[
            Column("site_id", "string", nullable=False),
            Column("forest_cover", "float", required=False),
            Column("humidity", "float", required=False),
            Column("temperature", "float", required=False),
            Column("precipitation", "float", required=False),
            Column("diet", "string", required=False),
        ],
    ),
    "crosswalk": TableSchema(
        name="crosswalk",
        unit="identifier relationship",
        primary_key=("entity", "source_id", "canonical_id"),
        columns=[
            Column("entity", "string", nullable=False),
            Column("source_id", "string", nullable=False),
            Column("canonical_id", "string", nullable=False),
            Column("source_table", "string", required=False),
            Column("source_column", "string", required=False),
            Column("evidence", "string", required=False),
            Column("provenance", "string", required=False),
        ],
    ),
}


def _dtype_ok(series: pd.Series, dtype: str) -> bool:
    if dtype in {"string", "category"}:
        return series.dtype == object or isinstance(series.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(series)
    if dtype == "float":
        return pd.api.types.is_numeric_dtype(series)
    if dtype == "int":
        return pd.api.types.is_integer_dtype(series) or pd.api.types.is_float_dtype(series)
    if dtype == "bool":
        return pd.api.types.is_bool_dtype(series) or series.dropna().isin([True, False]).all()
    return True


def validate_table(df: pd.DataFrame, schema_name: str, strict: bool = True) -> List[str]:
    """Validate ``df`` against a canonical schema.

    Returns a list of problem descriptions. With ``strict=True`` a non-empty
    list raises :class:`SchemaError` instead of being returned.
    """
    schema = CANONICAL_SCHEMAS.get(schema_name)
    if schema is None:
        raise KeyError(f"unknown canonical schema {schema_name!r}")

    problems: List[str] = []

    for col in schema.columns:
        if col.name not in df.columns:
            if col.required:
                problems.append(f"{schema_name}: missing required column '{col.name}'")
            continue
        series = df[col.name]
        if not col.nullable and series.isna().any():
            n = int(series.isna().sum())
            problems.append(f"{schema_name}: column '{col.name}' has {n} null values but is not nullable")
        if not _dtype_ok(series, col.dtype):
            problems.append(
                f"{schema_name}: column '{col.name}' has dtype {series.dtype}, expected {col.dtype}"
            )
        if col.allowed is not None:
            observed = set(series.dropna().astype(str).unique())
            unexpected = sorted(observed - set(col.allowed))
            if unexpected:
                problems.append(
                    f"{schema_name}: column '{col.name}' has values outside "
                    f"{list(col.allowed)}: {unexpected[:8]}"
                )

    missing_pk = [k for k in schema.primary_key if k not in df.columns]
    if missing_pk:
        problems.append(f"{schema_name}: primary key columns absent: {missing_pk}")
    elif len(df):
        dup = df.duplicated(subset=list(schema.primary_key)).sum()
        if dup:
            problems.append(
                f"{schema_name}: primary key {list(schema.primary_key)} is not unique "
                f"({int(dup)} duplicate rows)"
            )

    if "detection" in df.columns and "y" in df.columns:
        mapping = {"positive": 1.0, "negative": 0.0}
        expected = df["detection"].astype(str).map(mapping)
        actual = pd.to_numeric(df["y"], errors="coerce")
        mismatch = int(((expected.notna()) & (expected != actual)).sum())
        if mismatch:
            problems.append(
                f"{schema_name}: {mismatch} rows where numeric 'y' disagrees with 'detection'"
            )
        leaked = int(((df["detection"].astype(str) == "unassessable") & actual.notna()).sum())
        if leaked:
            problems.append(
                f"{schema_name}: {leaked} unassessable rows carry a numeric response; "
                "unassessable cells must stay missing, not become zeros"
            )

    if strict and problems:
        raise SchemaError("; ".join(problems))
    return problems


def check_foreign_keys(tables: Dict[str, pd.DataFrame], strict: bool = False) -> List[str]:
    """Check declared foreign keys across the harmonized tables."""
    problems: List[str] = []
    for name, schema in CANONICAL_SCHEMAS.items():
        child = tables.get(name)
        if child is None:
            continue
        for column, target in schema.foreign_keys.items():
            parent_name, parent_col = target.split(".")
            parent = tables.get(parent_name)
            if parent is None or column not in child.columns or parent_col not in parent.columns:
                continue
            known = set(parent[parent_col].dropna().astype(str))
            values = child[column].dropna().astype(str)
            orphans = sorted(set(values) - known)
            if orphans:
                problems.append(
                    f"{name}.{column}: {len(orphans)} values absent from {target} "
                    f"(e.g. {orphans[:5]})"
                )
    if strict and problems:
        raise SchemaError("; ".join(problems))
    return problems


def schema_as_dict() -> Dict[str, Any]:
    """Machine-readable dump used to generate ``docs/DATA_DICTIONARY.md``."""
    return {
        name: {
            "unit": s.unit,
            "primary_key": list(s.primary_key),
            "foreign_keys": dict(s.foreign_keys),
            "columns": [
                {
                    "name": c.name,
                    "dtype": c.dtype,
                    "required": c.required,
                    "nullable": c.nullable,
                    "allowed": list(c.allowed) if c.allowed else None,
                    "description": c.description,
                }
                for c in s.columns
            ],
        }
        for name, s in CANONICAL_SCHEMAS.items()
    }
