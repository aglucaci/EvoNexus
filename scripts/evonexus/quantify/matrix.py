"""Turn the canonical tables into one analysis frame of eligible cells.

Eligibility is applied once, here, so that every model in the M0-M5 ladder and
every sensitivity variant sees the same observation set and their predictive
comparisons stay well defined.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from evonexus.core.config import load_config
from evonexus.core.io import write_json, write_table
from evonexus.core.provenance import RunRecord

ECOLOGY_COLUMNS = ("forest_cover", "humidity", "temperature", "precipitation")


def _standardize(values: pd.Series) -> tuple[pd.Series, float, float]:
    numeric = pd.to_numeric(values, errors="coerce")
    mu = float(np.nanmean(numeric)) if numeric.notna().any() else 0.0
    sd = float(np.nanstd(numeric)) if numeric.notna().any() else 0.0
    sd = sd if sd > 0 else 1.0
    return (numeric - mu) / sd, mu, sd


def build_analysis_frame(
    tables: dict[str, pd.DataFrame],
    filters: dict[str, Any],
    resolution: str = "lineage",
    drop_mixed_pools: bool = True,
    require_known_habitat: bool = True,
) -> dict[str, Any]:
    """Join the canonical tables into one eligible-cell frame.

    Returns the frame, the standardization constants (needed to interpret and
    to re-apply at prediction time) and an exclusion ledger.
    """
    lib, hosts, sites, vir, det = (
        tables["libraries"].copy(), tables["hosts"].copy(), tables["sites"].copy(),
        tables["viruses"].copy(), tables["detections"].copy(),
    )
    eco = tables.get("ecology", pd.DataFrame()).copy()
    ledger: list[dict[str, Any]] = []

    def note(stage: str, n_cells: int, reason: str) -> None:
        ledger.append({"stage": stage, "cells_remaining": int(n_cells), "reason": reason})

    if det.empty or lib.empty or vir.empty:
        return {"frame": pd.DataFrame(), "standardization": {}, "ledger": ledger,
                "status": "empty_input"}

    note("start", len(det), "All harmonized library x lineage cells")

    # --- library-level attributes --------------------------------------
    lib["library_id"] = lib["library_id"].astype(str)
    lib = lib.set_index("library_id")
    site_habitat = sites.set_index("site_id")["habitat"] if len(sites) else pd.Series(dtype=object)
    host_species = (
        hosts.set_index("host_id")["host_species_resolved"] if len(hosts) else pd.Series(dtype=object)
    )
    host_mixed = (
        hosts.set_index("host_id")["mixed_species_pool"] if len(hosts) else pd.Series(dtype=bool)
    )
    host_poolsize = hosts.set_index("host_id")["pool_size"] if len(hosts) else pd.Series(dtype=float)

    lib_attr = pd.DataFrame(index=lib.index)
    lib_attr["site_id"] = lib["site_id"].astype(str)
    lib_attr["habitat"] = lib_attr["site_id"].map(site_habitat).fillna("unknown")
    lib_attr["host_id"] = lib["host_id"].astype(str)
    lib_attr["host_species"] = lib_attr["host_id"].map(host_species).fillna("unknown")
    lib_attr["mixed_pool"] = lib_attr["host_id"].map(host_mixed).fillna(False).astype(bool)
    lib_attr["pool_size"] = lib_attr["host_id"].map(host_poolsize)
    lib_attr["tissue_group"] = lib["tissue_group"].fillna("unknown")
    lib_attr["sequencing_batch"] = lib["sequencing_batch"].fillna("unknown").astype(str)
    lib_attr["effective_depth"] = pd.to_numeric(lib["effective_depth"], errors="coerce")
    lib_attr["collection_time"] = pd.to_numeric(lib["collection_time_numeric"], errors="coerce")
    lib_attr["is_negative_control"] = lib["is_negative_control"].fillna(False).astype(bool)

    if len(eco):
        eco = eco.set_index("site_id")
        for col in ECOLOGY_COLUMNS:
            if col in eco.columns:
                lib_attr[col] = lib_attr["site_id"].map(eco[col])

    # --- virus-level attributes ----------------------------------------
    vir["virus_id"] = vir["virus_id"].astype(str)
    vir_attr = vir.set_index("virus_id")[["virus_family", "eligible"]]

    df = det.copy()
    df["library_id"] = df["library_id"].astype(str)
    df["virus_id"] = df["virus_id"].astype(str)
    df = df.join(lib_attr, on="library_id").join(vir_attr, on="virus_id")

    # --- eligibility filters -------------------------------------------
    df = df[~df["is_negative_control"].fillna(False)]
    note("drop_negative_controls", len(df), "Negative/blank control libraries removed")

    df = df[df["detection"].isin(["positive", "negative"])]
    note("assessable_only", len(df),
         "Unassessable cells dropped: no response, never imputed as zeros")

    df = df[df["eligible"].fillna(False)]
    note("eligible_lineages", len(df), "Mammal-associated RNA viral families only")

    if require_known_habitat:
        df = df[df["habitat"].isin(["settlement", "wild"])]
        note("known_habitat", len(df), "Libraries from sites with a classified habitat")

    if drop_mixed_pools:
        df = df[(~df["mixed_pool"]) & (~df["host_species"].isin(["unknown", "mixed"]))]
        note("single_host_species", len(df),
             "Mixed-species and unidentified pools excluded from the single-host model")

    df = df[df["effective_depth"].notna() & (df["effective_depth"] > 0)]
    note("depth_available", len(df), "Libraries with a usable effective-depth denominator")

    if df.empty:
        return {"frame": df, "standardization": {}, "ledger": ledger, "status": "empty_after_filters"}

    # --- lineage resolution --------------------------------------------
    rules = filters.get("lineage_eligibility", {}) or {}
    min_libs = int(rules.get("min_libraries_detected", 2))
    pool_rare = bool(rules.get("pool_rare_into_family", True))

    if resolution == "family":
        df["unit_id"] = df["virus_family"].astype(str)
        df = (
            df.groupby(["library_id", "unit_id"], as_index=False)
            .agg(
                y=("y", "max"),
                virus_family=("virus_family", "first"),
                **{c: (c, "first") for c in (
                    "habitat", "host_id", "host_species", "pool_size", "tissue_group",
                    "sequencing_batch", "site_id", "effective_depth", "collection_time",
                ) if c in df.columns},
                **{c: (c, "first") for c in ECOLOGY_COLUMNS if c in df.columns},
            )
        )
    else:
        if resolution == "strict":
            strict_mask = (df["detection"] == "positive") & (
                pd.to_numeric(df.get("coverage_breadth"), errors="coerce").fillna(0) < 0.5
            )
            df.loc[strict_mask, "y"] = 0.0
            df.loc[strict_mask, "detection"] = "negative"
            note("strict_resolution", len(df),
                 "Positives below 0.5 coverage breadth downgraded to negatives")
        pos_per_lineage = df[df["y"] > 0].groupby("virus_id")["library_id"].nunique()
        rare = set(pos_per_lineage[pos_per_lineage < min_libs].index) | (
            set(df["virus_id"].unique()) - set(pos_per_lineage.index)
        )
        if pool_rare:
            df["unit_id"] = np.where(
                df["virus_id"].isin(rare),
                "FAMILYPOOL:" + df["virus_family"].astype(str),
                df["virus_id"].astype(str),
            )
            note("pool_rare_lineages", len(df),
                 f"{len(rare)} lineages detected in <{min_libs} libraries pooled into their family")
        else:
            df = df[~df["virus_id"].isin(rare)]
            df["unit_id"] = df["virus_id"].astype(str)
            note("drop_rare_lineages", len(df),
                 f"{len(rare)} lineages detected in <{min_libs} libraries dropped")
        df = df.groupby(["library_id", "unit_id"], as_index=False).agg(
            y=("y", "max"),
            virus_family=("virus_family", "first"),
            **{c: (c, "first") for c in (
                "habitat", "host_id", "host_species", "pool_size", "tissue_group",
                "sequencing_batch", "site_id", "effective_depth", "collection_time",
            ) if c in df.columns},
            **{c: (c, "first") for c in ECOLOGY_COLUMNS if c in df.columns},
        )

    # --- model covariates ------------------------------------------------
    df["U"] = (df["habitat"] == "settlement").astype(float)
    df["log_depth"] = np.log(pd.to_numeric(df["effective_depth"], errors="coerce"))
    df["log_pool_size"] = np.log(pd.to_numeric(df["pool_size"], errors="coerce").fillna(1.0).clip(lower=1.0))

    standardization: dict[str, dict[str, float]] = {}
    for col, target in (("log_depth", "log_depth_z"), ("log_pool_size", "log_pool_size_z"),
                        ("collection_time", "time_z")):
        values, mu, sd = _standardize(df[col])
        df[target] = values.fillna(0.0)
        standardization[target] = {"source": col, "mean": mu, "sd": sd}

    ecology_available = []
    for col in ECOLOGY_COLUMNS:
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().mean() > 0.8:
            values, mu, sd = _standardize(df[col])
            df[f"{col}_z"] = values.fillna(0.0)
            standardization[f"{col}_z"] = {"source": col, "mean": mu, "sd": sd}
            ecology_available.append(f"{col}_z")

    df = df.rename(columns={"unit_id": "virus_unit"})
    keep = [
        "library_id", "virus_unit", "virus_family", "y", "habitat", "U",
        "host_species", "host_id", "tissue_group", "site_id", "sequencing_batch",
        "effective_depth", "log_depth_z", "log_pool_size_z", "time_z",
    ] + ecology_available
    frame = df[[c for c in keep if c in df.columns]].reset_index(drop=True)
    frame["y"] = pd.to_numeric(frame["y"], errors="coerce")

    note("final", len(frame), "Cells entering the likelihood")

    summary = {
        "n_cells": len(frame),
        "n_libraries": int(frame["library_id"].nunique()),
        "n_virus_units": int(frame["virus_unit"].nunique()),
        "n_families": int(frame["virus_family"].nunique()),
        "n_host_species": int(frame["host_species"].nunique()),
        "n_sites": int(frame["site_id"].nunique()),
        "n_batches": int(frame["sequencing_batch"].nunique()),
        "detection_rate": float(frame["y"].mean()),
        "detection_rate_by_habitat": frame.groupby("habitat")["y"].mean().to_dict(),
        "libraries_by_habitat": frame.groupby("habitat")["library_id"].nunique().to_dict(),
        "ecology_covariates": ecology_available,
        "resolution": resolution,
    }

    return {
        "frame": frame,
        "standardization": standardization,
        "ledger": ledger,
        "summary": summary,
        "ecology_available": ecology_available,
        "status": "ok",
    }


# ---------------------------------------------------------------------------
# sensitivity restrictions
# ---------------------------------------------------------------------------
def _restrict_overlap_only(frame: pd.DataFrame) -> pd.DataFrame:
    table = pd.crosstab(frame["host_species"], frame["habitat"])
    if not {"settlement", "wild"}.issubset(table.columns):
        return frame.iloc[0:0]
    keep = table.index[(table["settlement"] > 0) & (table["wild"] > 0)]
    return frame[frame["host_species"].isin(keep)]


def _restrict_no_pools(frame: pd.DataFrame, tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    hosts = tables.get("hosts", pd.DataFrame())
    if hosts.empty:
        return frame
    singles = set(hosts.loc[hosts["pool_status"].astype(str) == "individual", "host_id"].astype(str))
    if "host_id" not in frame.columns:
        return frame
    return frame[frame["host_id"].astype(str).isin(singles)]


def _restrict_drop_rare(frame: pd.DataFrame, min_detections: int = 5) -> pd.DataFrame:
    counts = frame[frame["y"] > 0].groupby("virus_unit").size()
    keep = set(counts[counts >= min_detections].index)
    return frame[frame["virus_unit"].isin(keep)]


def _restrict_depth_matched(frame: pd.DataFrame) -> pd.DataFrame:
    if "habitat" not in frame or frame["habitat"].nunique() < 2:
        return frame
    bounds = frame.groupby("habitat")["log_depth_z"].agg(["min", "max"])
    lo = float(bounds["min"].max())
    hi = float(bounds["max"].min())
    if lo >= hi:
        return frame.iloc[0:0]
    return frame[(frame["log_depth_z"] >= lo) & (frame["log_depth_z"] <= hi)]


RESTRICTIONS: dict[str, Callable[..., pd.DataFrame]] = {
    "none": lambda frame, tables=None: frame,
    "overlap_only": lambda frame, tables=None: _restrict_overlap_only(frame),
    "no_pools": lambda frame, tables=None: _restrict_no_pools(frame, tables or {}),
    "drop_rare_lineages": lambda frame, tables=None: _restrict_drop_rare(frame),
    "depth_matched": lambda frame, tables=None: _restrict_depth_matched(frame),
}


def apply_restriction(frame: pd.DataFrame, name: str,
                      tables: dict[str, pd.DataFrame] | None = None) -> pd.DataFrame:
    func = RESTRICTIONS.get(name)
    if func is None:
        raise KeyError(f"unknown restriction {name!r}; available: {sorted(RESTRICTIONS)}")
    return func(frame, tables).reset_index(drop=True)


def wide_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """Binary library x lineage matrix; unassessable cells stay NaN."""
    if frame.empty:
        return pd.DataFrame()
    return (
        frame.pivot_table(index="library_id", columns="virus_unit", values="y", aggfunc="max")
        .reset_index()
    )


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def run_matrix(
    processed_dir: str | Path,
    out_dir: str | Path,
    resolution: str = "lineage",
    filters_config: str = "viral_filters",
    synthetic: bool = False,
) -> dict[str, Any]:
    from evonexus.qc.run import load_canonical

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    record = RunRecord(step="matrix", params={"resolution": resolution, "synthetic": synthetic})
    record.synthetic = synthetic

    tables = load_canonical(processed_dir)
    filters = load_config(filters_config)
    built = build_analysis_frame(tables, filters, resolution=resolution)
    frame = built["frame"]

    frame_path = write_table(frame, out / f"analysis_frame_{resolution}")
    wide_path = write_table(wide_matrix(frame), out / f"detection_matrix_{resolution}")
    ledger_path = write_table(pd.DataFrame(built["ledger"]), out / f"matrix_ledger_{resolution}")
    meta_path = write_json(
        {
            "status": built["status"],
            "summary": built.get("summary", {}),
            "standardization": built.get("standardization", {}),
            "ecology_available": built.get("ecology_available", []),
            "synthetic": synthetic,
        },
        out / f"matrix_metadata_{resolution}.json",
    )
    for p in (frame_path, wide_path, ledger_path, meta_path):
        record.add_output(p, role="matrix")

    record.note(
        f"resolution={resolution}: {len(frame)} cells, "
        f"{frame['library_id'].nunique() if len(frame) else 0} libraries, "
        f"{frame['virus_unit'].nunique() if len(frame) else 0} lineage units"
    )
    record.finish(out / f"run_record_matrix_{resolution}.json",
                  status="ok" if built["status"] == "ok" else "warning")

    return {
        "frame": str(frame_path),
        "matrix": str(wide_path),
        "ledger": str(ledger_path),
        "metadata": str(meta_path),
        "status": built["status"],
        "summary": built.get("summary", {}),
    }
