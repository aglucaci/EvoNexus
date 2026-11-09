"""Harmonization and schema tests.

These guard the claims the rest of the pipeline depends on: that unknown is
not zero, that a pooled library is not an animal, and that a mapping is never
applied without being resolved.
"""

from __future__ import annotations

import pandas as pd
import pytest

from evonexus.core.config import load_config, resolve_variant
from evonexus.core.schema import (
    CANONICAL_SCHEMAS,
    SchemaError,
    check_foreign_keys,
    validate_table,
)
from evonexus.harmonize.mapping import resolve_mapping


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------
def test_every_config_loads():
    for name in ("datasets", "metadata_mapping", "viral_filters", "priors",
                 "model_comparisons", "figures", "qc_gates", "evolution"):
        assert load_config(name), f"{name}.yaml is empty or missing"


def test_real_mapping_is_not_preapproved():
    """A shipped mapping must never claim human approval it did not get."""
    assert load_config("metadata_mapping").get("approved") is False


def test_prior_variants_resolve_and_inherit():
    priors = load_config("priors")
    base = resolve_variant(priors, "regularizing")
    tight = resolve_variant(priors, "tight")
    assert base["habitat_beta"]["sigma"] == 1.0
    assert tight["habitat_beta"]["sigma"] == 0.5
    # Inherited keys survive the override.
    assert tight["depth_gamma"] == base["depth_gamma"]
    assert tight["_inheritance_chain"] == ["regularizing", "tight"]


def test_unknown_variant_raises():
    with pytest.raises(KeyError):
        resolve_variant(load_config("priors"), "nonexistent")


# ---------------------------------------------------------------------------
# mapping resolution
# ---------------------------------------------------------------------------
def test_unapproved_mapping_is_blocked(inventory):
    mapping = resolve_mapping(inventory, load_config("metadata_mapping"),
                              allow_role_guess=False, require_approval=True)
    assert mapping.approved is False
    assert mapping.warnings
    assert mapping.blocking_problems, "an unresolved mapping must block harmonization"


def test_mapping_reports_every_field(inventory, synthetic_raw):
    mapping = resolve_mapping(
        inventory, load_config("metadata_mapping"),
        overrides_path=synthetic_raw / "synthetic" / "mapping_overrides.json",
        allow_role_guess=True, require_approval=False,
    )
    libraries = mapping.entities["libraries"]
    assert libraries.source_path is not None
    resolved = {f.canonical: f.source_column for f in libraries.fields}
    assert resolved["library_id"] == "Library"
    assert resolved["read_depth"] == "clean_reads"
    assert not libraries.missing_required


# ---------------------------------------------------------------------------
# canonical tables
# ---------------------------------------------------------------------------
def test_all_canonical_tables_validate(tables):
    for name in CANONICAL_SCHEMAS:
        assert name in tables, f"{name} was not produced"
        validate_table(tables[name], name, strict=True)


def test_foreign_keys_hold(tables):
    assert check_foreign_keys(tables, strict=False) == []


def test_habitat_has_an_explicit_unknown_level(tables):
    sites = tables["sites"]
    assert set(sites["habitat"]).issubset({"settlement", "wild", "unknown"})
    # The simulator emits blank habitats; they must land in `unknown`.
    assert (sites["habitat"] == "unknown").any()


def test_mixed_species_pools_are_flagged_not_split(tables):
    hosts = tables["hosts"]
    mixed = hosts[hosts["mixed_species_pool"].fillna(False)]
    assert len(mixed) > 0
    assert (mixed["host_species_resolved"] == "mixed").all()
    # One row per pool: no duplication into one row per constituent species.
    assert hosts["host_id"].is_unique


def test_unassessable_cells_carry_no_response(tables):
    det = tables["detections"]
    unassessable = det[det["detection"] == "unassessable"]
    assert len(unassessable) > 0, "the fixture should contain unscreened cells"
    assert unassessable["y"].isna().all(), "unassessable cells must not become zeros"


def test_positives_from_unscreened_libraries_are_excluded_but_recorded(tables):
    """The asymmetry that would bias detection upward at low depth."""
    det = tables["detections"]
    outside = det[det["evidence"] == "positive_outside_assessed_universe"]
    if len(outside):
        assert (outside["detection"] == "unassessable").all()
        assert outside["y"].isna().all()
        assert (outside["flag"].astype(str) != "").all()


def test_detection_counts_are_consistent(tables):
    det = tables["detections"]
    assert (det.loc[det["detection"] == "positive", "y"] == 1.0).all()
    assert (det.loc[det["detection"] == "negative", "y"] == 0.0).all()
    assert det["detection"].isin(["positive", "negative", "unassessable"]).all()


def test_ineligible_viral_families_are_excluded(tables):
    viruses = tables["viruses"]
    insect = viruses[viruses["virus_family"] == "Dicistroviridae"]
    if len(insect):
        assert not insect["eligible"].any()
        assert insect["eligibility_reason"].str.contains("excluded").all()


def test_crosswalk_records_source_provenance(tables):
    cw = tables["crosswalk"]
    assert len(cw) > 0
    assert set(["entity", "source_id", "canonical_id", "source_table"]).issubset(cw.columns)
    assert cw["source_table"].notna().any()


def test_canonical_tables_survive_concat(tables):
    """Regression: a DataFrame hidden in `.attrs` broke pandas concat.

    `pd.concat` compares the `attrs` of its inputs, and a DataFrame stored
    there raises on its ambiguous truth value. The source frame is now passed
    explicitly instead. This only surfaced with pyarrow installed, so it was
    invisible in a CSV-fallback environment.
    """
    for name, df in tables.items():
        if not len(df):
            continue
        assert df.attrs == {}, f"{name} carries unexpected DataFrame.attrs"
        doubled = pd.concat([df, df], ignore_index=True)
        assert len(doubled) == 2 * len(df)


def test_effective_depth_definition_is_uniform(tables):
    libs = tables["libraries"]
    assert libs["effective_depth_definition"].nunique() == 1


# ---------------------------------------------------------------------------
# schema enforcement
# ---------------------------------------------------------------------------
def test_schema_rejects_zeros_smuggled_into_unassessable_cells():
    df = pd.DataFrame(
        {
            "library_id": ["L1"],
            "virus_id": ["V1"],
            "detection": ["unassessable"],
            "y": [0.0],
        }
    )
    with pytest.raises(SchemaError, match="unassessable"):
        validate_table(df, "detections", strict=True)


def test_schema_rejects_duplicate_primary_keys():
    df = pd.DataFrame({"site_id": ["S1", "S1"], "habitat": ["wild", "wild"]})
    problems = validate_table(df, "sites", strict=False)
    assert any("not unique" in p for p in problems)


def test_schema_rejects_unknown_habitat_level():
    df = pd.DataFrame({"site_id": ["S1"], "habitat": ["suburban"]})
    problems = validate_table(df, "sites", strict=False)
    assert any("outside" in p for p in problems)
