"""Posterior standardization and the primary habitat estimand.

The estimand
------------
The difference in *average lineage-level detection probability* between
settlement and wild habitats, standardized to one common sampling
distribution: the same host species, tissues, sequencing depths, pool sizes
and collection times, at a *new* site and a *new* sequencing batch drawn from
their estimated population distributions.

What it is not
--------------
It is not a causal effect of urbanization. Habitat is not randomized, sites
are not interchangeable, and the standardization removes only the measured
and modelled differences in sampling. Every artifact this module writes
carries that statement.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np
import pandas as pd

from evonexus.bayesian.diagnostics import hdi
from evonexus.bayesian.predict import PosteriorPredictor, expit

ESTIMAND_CAVEAT = (
    "Observational association standardized to a common sampling distribution. "
    "Not a causal effect of urbanization."
)

_PROFILE_COLUMNS = (
    "host_species", "tissue_group", "log_depth_z", "log_pool_size_z", "time_z",
)


def reference_population(
    frame: pd.DataFrame,
    overlap_only: bool = True,
    max_libraries: int = 120,
    seed: int = 20250101,
) -> pd.DataFrame:
    """The library profiles that the contrast is standardized over.

    Restricting to host species observed in both habitats is deliberate: it
    keeps the contrast anchored on a comparison that the data can actually
    support, rather than extrapolating a species seen only in one habitat.
    """
    if frame.empty:
        return frame
    libs = frame.drop_duplicates(subset=["library_id"]).copy()
    if overlap_only:
        table = pd.crosstab(libs["host_species"], libs["habitat"])
        if {"settlement", "wild"}.issubset(table.columns):
            keep = table.index[(table["settlement"] > 0) & (table["wild"] > 0)]
            if len(keep):
                libs = libs[libs["host_species"].isin(keep)]
    if len(libs) > max_libraries:
        libs = libs.sample(max_libraries, random_state=seed)
    ecology = [c for c in frame.columns if c.endswith("_z")
               and c not in {"log_depth_z", "log_pool_size_z", "time_z"}]
    columns = [c for c in _PROFILE_COLUMNS if c in libs.columns] + ecology + ["library_id", "habitat"]
    return libs[columns].reset_index(drop=True)


def _lineage_table(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame[["virus_unit", "virus_family"]]
        .drop_duplicates()
        .sort_values("virus_unit")
        .reset_index(drop=True)
    )


def _build_grid(profiles: pd.DataFrame, lineages: pd.DataFrame) -> pd.DataFrame:
    grid = profiles.merge(lineages, how="cross")
    grid["site_id"] = "__new_site__"
    grid["sequencing_batch"] = "__new_batch__"
    return grid


def standardized_predictions(
    posterior: Dict[str, Any],
    frame: pd.DataFrame,
    overlap_only: bool = True,
    max_libraries: int = 120,
    seed: int = 20250101,
    chunk_size: int = 4000,
) -> Dict[str, Any]:
    """Posterior draws of standardized detection probability per habitat."""
    predictor = PosteriorPredictor(posterior)
    profiles = reference_population(frame, overlap_only, max_libraries, seed)
    lineages = _lineage_table(frame)
    if profiles.empty or lineages.empty:
        return {"status": "empty", "reason": "no reference population or no lineages"}

    grid = _build_grid(profiles, lineages)
    families = grid["virus_family"].astype(str).to_numpy()
    unique_families = sorted(set(families))
    hosts = grid["host_species"].astype(str).to_numpy()
    unique_hosts = sorted(set(hosts))

    S = predictor.n_samples
    totals = {u: np.zeros(S) for u in (0.0, 1.0)}
    family_totals = {u: {f: np.zeros(S) for f in unique_families} for u in (0.0, 1.0)}
    host_totals = {u: {h: np.zeros(S) for h in unique_hosts} for u in (0.0, 1.0)}
    family_counts = {f: int((families == f).sum()) for f in unique_families}
    host_counts = {h: int((hosts == h).sum()) for h in unique_hosts}
    n_cells = len(grid)

    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        block = grid.iloc[start:stop]
        fam_chunk = families[start:stop]
        host_chunk = hosts[start:stop]
        for u in (0.0, 1.0):
            # One RNG per (chunk, habitat) so the same new-site draws are used
            # for both habitats: the contrast is then paired, not independent.
            rng = np.random.default_rng(seed + start)
            p = expit(predictor.eta(block, integrate=("site", "batch"), rng=rng,
                                    habitat_override=u))
            totals[u] += p.sum(axis=1)
            for f in unique_families:
                mask = fam_chunk == f
                if mask.any():
                    family_totals[u][f] += p[:, mask].sum(axis=1)
            for h in unique_hosts:
                mask = host_chunk == h
                if mask.any():
                    host_totals[u][h] += p[:, mask].sum(axis=1)

    p_settlement = totals[1.0] / n_cells
    p_wild = totals[0.0] / n_cells

    return {
        "status": "ok",
        "p_settlement": p_settlement,
        "p_wild": p_wild,
        "difference": p_settlement - p_wild,
        "risk_ratio": np.where(p_wild > 0, p_settlement / np.maximum(p_wild, 1e-12), np.nan),
        "by_family": {
            f: {
                "p_settlement": family_totals[1.0][f] / max(family_counts[f], 1),
                "p_wild": family_totals[0.0][f] / max(family_counts[f], 1),
            }
            for f in unique_families
        },
        "by_host": {
            h: {
                "p_settlement": host_totals[1.0][h] / max(host_counts[h], 1),
                "p_wild": host_totals[0.0][h] / max(host_counts[h], 1),
            }
            for h in unique_hosts
        },
        "reference": {
            "n_library_profiles": int(len(profiles)),
            "n_lineages": int(len(lineages)),
            "n_standardization_cells": n_cells,
            "overlap_only": overlap_only,
            "libraries_used": profiles["library_id"].tolist(),
            "host_species_used": sorted(set(profiles["host_species"].astype(str))),
        },
        "engine": predictor.engine,
        "caveat": ESTIMAND_CAVEAT,
    }


def _summary_row(label: str, draws: np.ndarray, scale: str, hdi_prob: float = 0.95) -> Dict[str, Any]:
    flat = np.asarray(draws, dtype=float).reshape(-1)
    lo, hi = hdi(flat, hdi_prob)
    return {
        "quantity": label,
        "scale": scale,
        "mean": float(np.mean(flat)),
        "median": float(np.median(flat)),
        "sd": float(np.std(flat, ddof=1)) if flat.size > 1 else float("nan"),
        f"hdi_{int(hdi_prob * 100)}_low": lo,
        f"hdi_{int(hdi_prob * 100)}_high": hi,
        "p_direction_positive": float(np.mean(flat > 0)),
        "n_draws": int(flat.size),
    }


def habitat_contrast(
    posterior: Dict[str, Any],
    frame: pd.DataFrame,
    hdi_prob: float = 0.95,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Standardized settlement-versus-wild contrast, overall and per family."""
    std = standardized_predictions(posterior, frame, **kwargs)
    if std.get("status") != "ok":
        return {"status": std.get("status", "failed"), "reason": std.get("reason"),
                "table": pd.DataFrame(), "caveat": ESTIMAND_CAVEAT}

    rows: List[Dict[str, Any]] = [
        _summary_row("standardized detection probability (settlement)",
                     std["p_settlement"], "probability", hdi_prob),
        _summary_row("standardized detection probability (wild)",
                     std["p_wild"], "probability", hdi_prob),
        _summary_row("standardized difference (settlement - wild)",
                     std["difference"], "probability difference", hdi_prob),
    ]

    predictor = PosteriorPredictor(posterior)
    beta_u = predictor.fixed("habitat")
    if np.any(beta_u != 0):
        rows.append(_summary_row("habitat coefficient beta_U", beta_u, "log odds", hdi_prob))
        rows.append(_summary_row("habitat odds ratio", np.exp(beta_u), "odds ratio", hdi_prob))

    family_rows: List[Dict[str, Any]] = []
    for family, draws in std["by_family"].items():
        diff = draws["p_settlement"] - draws["p_wild"]
        row = _summary_row(f"difference — {family}", diff, "probability difference", hdi_prob)
        row["virus_family"] = family
        family_rows.append(row)

    host_rows: List[Dict[str, Any]] = []
    for host, draws in std["by_host"].items():
        diff = draws["p_settlement"] - draws["p_wild"]
        row = _summary_row(f"difference — {host}", diff, "probability difference", hdi_prob)
        row["host_species"] = host
        host_rows.append(row)

    return {
        "status": "ok",
        "table": pd.DataFrame(rows),
        "by_family": pd.DataFrame(family_rows),
        "by_host": pd.DataFrame(host_rows),
        "draws": {
            "p_settlement": std["p_settlement"],
            "p_wild": std["p_wild"],
            "difference": std["difference"],
        },
        "reference": std["reference"],
        "engine": std["engine"],
        "caveat": ESTIMAND_CAVEAT,
    }


def raw_detection_rates(frame: pd.DataFrame) -> pd.DataFrame:
    """Unadjusted per-library detection rates, for the comparison panel."""
    if frame.empty:
        return pd.DataFrame()
    per_library = (
        frame.groupby(["library_id", "habitat", "host_species"], as_index=False)
        .agg(
            detection_rate=("y", "mean"),
            n_detected=("y", "sum"),
            n_lineages_tested=("y", "size"),
        )
    )
    return per_library


def posterior_predictive_replicates(
    posterior: Dict[str, Any],
    frame: pd.DataFrame,
    n_replicates: int = 200,
    seed: int = 20250101,
) -> Dict[str, Any]:
    """Replicate datasets under the fitted model, for posterior checks."""
    predictor = PosteriorPredictor(posterior)
    if frame.empty:
        return {"status": "empty"}
    rng = np.random.default_rng(seed)
    idx = rng.choice(predictor.n_samples, size=min(n_replicates, predictor.n_samples),
                     replace=False)
    p = expit(predictor.eta(frame, integrate=(), rng=np.random.default_rng(seed)))[idx]
    y_rep = (rng.random(p.shape) < p).astype(float)
    observed = pd.to_numeric(frame["y"], errors="coerce").to_numpy(float)

    def summarize(matrix: np.ndarray, column: str) -> Dict[str, Any]:
        groups = frame[column].astype(str).to_numpy()
        levels = sorted(set(groups))
        return {
            lvl: matrix[:, groups == lvl].mean(axis=1).tolist()
            for lvl in levels
        }

    return {
        "status": "ok",
        "observed_detection_rate": float(observed.mean()),
        "replicate_detection_rate": y_rep.mean(axis=1).tolist(),
        "observed_by_habitat": frame.groupby("habitat")["y"].mean().to_dict(),
        "replicate_by_habitat": summarize(y_rep, "habitat"),
        "observed_lineages_per_library": float(
            frame.groupby("library_id")["y"].sum().mean()
        ),
        "replicate_lineages_per_library": [
            float(pd.Series(row).groupby(frame["library_id"].to_numpy()).sum().mean())
            for row in y_rep[: min(50, y_rep.shape[0])]
        ],
        "p_fitted_mean": p.mean(axis=0).tolist(),
        "n_replicates": int(y_rep.shape[0]),
    }
