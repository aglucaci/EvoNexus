"""Manuscript tables, emitted as data (CSV/Parquet) and as formatted Markdown.

A table whose inputs are unavailable is not written. It is listed in the
report as "not generated", with the reason.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import numpy as np
import pandas as pd

from evonexus.core.io import write_table, write_text
from evonexus.reporting.artifacts import ArtifactStore

Result = Tuple[pd.DataFrame | None, str]


def _fmt(value: Any, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "—"
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    if isinstance(value, (float, np.floating)):
        return f"{value:.{digits}g}"
    return str(value)


def to_markdown(df: pd.DataFrame, digits: int = 3, max_rows: int = 60) -> str:
    if df is None or not len(df):
        return "_No rows._"
    shown = df.head(max_rows)
    header = "| " + " | ".join(str(c) for c in shown.columns) + " |"
    sep = "| " + " | ".join("---" for _ in shown.columns) + " |"
    lines = [header, sep]
    for _, row in shown.iterrows():
        lines.append("| " + " | ".join(_fmt(v, digits) for v in row) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_{len(df) - max_rows} further rows omitted; see the data file._")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# table builders
# ---------------------------------------------------------------------------
def table1(store: ArtifactStore) -> Result:
    libs = store.canonical("libraries")
    hosts = store.canonical("hosts")
    sites = store.canonical("sites")
    eco = store.canonical("ecology")
    if libs is None or not len(libs):
        return None, "canonical library table unavailable"

    merged = libs.merge(hosts[["host_id", "host_species_resolved", "pool_status", "pool_size"]],
                        on="host_id", how="left")
    merged = merged.merge(sites[["site_id", "habitat", "region"]], on="site_id", how="left")
    if eco is not None and len(eco):
        merged = merged.merge(eco, on="site_id", how="left")

    rows: List[Dict[str, Any]] = []
    for habitat, group in merged.groupby(merged["habitat"].fillna("unknown")):
        depth = pd.to_numeric(group["effective_depth"], errors="coerce")
        row = {
            "habitat": habitat,
            "n_libraries": int(len(group)),
            "n_sites": int(group["site_id"].nunique()),
            "n_source_pools": int(group["host_id"].nunique()),
            "n_host_species": int(
                group.loc[~group["host_species_resolved"].isin(["unknown", "mixed"]),
                          "host_species_resolved"].nunique()
            ),
            "pooled_libraries": int((group["pool_status"] == "pool").sum()),
            "median_effective_depth": float(depth.median()) if depth.notna().any() else np.nan,
            "iqr_effective_depth": (
                f"{depth.quantile(0.25):.3g}–{depth.quantile(0.75):.3g}"
                if depth.notna().any() else "—"
            ),
            "tissues": ", ".join(sorted(set(group["tissue_group"].dropna().astype(str)))[:5]),
        }
        for col in ("forest_cover", "humidity", "temperature", "precipitation"):
            if col in group:
                values = pd.to_numeric(group[col], errors="coerce")
                row[f"mean_{col}"] = float(values.mean()) if values.notna().any() else np.nan
        rows.append(row)

    total = {
        "habitat": "all",
        "n_libraries": int(len(merged)),
        "n_sites": int(merged["site_id"].nunique()),
        "n_source_pools": int(merged["host_id"].nunique()),
        "n_host_species": int(
            merged.loc[~merged["host_species_resolved"].isin(["unknown", "mixed"]),
                       "host_species_resolved"].nunique()
        ),
        "pooled_libraries": int((merged["pool_status"] == "pool").sum()),
    }
    rows.append(total)
    return pd.DataFrame(rows), "ok"


def table2(store: ArtifactStore) -> Result:
    contrast = store.contrast()
    if contrast is None or not len(contrast):
        return None, "habitat contrast unavailable (model not fitted or QC gate failed)"
    by_family = store.contrast_by_family()
    out = contrast.copy()
    out["stratum"] = "overall"
    if by_family is not None and len(by_family):
        fam = by_family.copy()
        fam["stratum"] = "viral family"
        out = pd.concat([out, fam], ignore_index=True)
    meta = store.contrast_metadata() or {}
    out["model"] = meta.get("model")
    out["engine"] = meta.get("engine")
    out["caveat"] = meta.get("caveat")
    return out, "ok"


def table3(store: ArtifactStore) -> Result:
    conn = store.connectivity()
    if conn is None or not len(conn):
        return None, "connectivity results unavailable"
    perm = store.connectivity_permutation()
    out = conn.copy()
    if perm is not None and len(perm):
        out = out.merge(perm, on=["host_i", "host_j"], how="left")
    meta = store.connectivity_metadata() or {}
    out["caveat"] = meta.get("caveat")
    return out, "ok"


def table4(store: ArtifactStore) -> Result:
    comparison = store.model_comparison()
    folds = store.cv_folds()
    if comparison is None and folds is None:
        return None, "neither model comparison nor cross-validation is available"
    frames: List[pd.DataFrame] = []
    if comparison is not None and len(comparison):
        comp = comparison.copy()
        comp["section"] = "in-sample model comparison (LOO)"
        frames.append(comp)
    if folds is not None and len(folds):
        agg = (
            folds.groupby("grouping")
            .agg(n_folds=("fold", "nunique"), mean_lpd=("mean_lpd", "mean"),
                 brier=("brier", "mean"), auc=("auc", "mean"), n=("n", "sum"))
            .reset_index()
        )
        agg["section"] = "held-out cross-validation"
        frames.append(agg)
    return pd.concat(frames, ignore_index=True), "ok"


def table_s1(store: ArtifactStore) -> Result:
    manifest = store.manifest()
    qc = store.qc()
    waterfall = store.waterfall()
    if manifest is None and qc is None:
        return None, "acquisition manifest and QC report are both unavailable"
    rows: List[Dict[str, Any]] = []
    for source in (manifest or {}).get("sources", []):
        rows.append(
            {
                "section": "data source",
                "item": source.get("source"),
                "detail": source.get("kind"),
                "status": source.get("status"),
                "value": source.get("accession") or (source.get("article_metadata") or {}).get("doi"),
            }
        )
    for gate in (qc or {}).get("gates", []):
        rows.append(
            {
                "section": "QC gate",
                "item": f"{gate['gate_id']} — {gate['name']}",
                "detail": gate.get("condition"),
                "status": gate.get("status"),
                "value": "; ".join(gate.get("failures") or []) or "conditions met",
            }
        )
    if waterfall is not None and len(waterfall):
        for _, r in waterfall.iterrows():
            rows.append(
                {
                    "section": "exclusion",
                    "item": r["stage"],
                    "detail": r["note"],
                    "status": "retained",
                    "value": f"{int(r['cells'])} cells, {int(r['libraries'])} libraries",
                }
            )
    return pd.DataFrame(rows), "ok"


def table_s2(store: ArtifactStore) -> Result:
    prior = store.prior_check()
    recovery = store.parameter_recovery()
    if prior is None and recovery is None:
        return None, "prior predictive check and parameter recovery are both unavailable"
    rows: List[Dict[str, Any]] = []
    if prior:
        for q, v in (prior.get("detection_rate_quantiles") or {}).items():
            rows.append({"section": "prior predictive", "item": f"detection rate {q}",
                         "value": v, "note": ""})
        rows.append({"section": "prior predictive", "item": "fraction in plausible band",
                     "value": prior.get("fraction_in_plausible_band"),
                     "note": "; ".join(prior.get("warnings") or []) or "no warnings"})
    if recovery is not None and len(recovery):
        for _, r in recovery.iterrows():
            rows.append(
                {
                    "section": "parameter recovery (synthetic)",
                    "item": r["parameter"],
                    "value": r["posterior_mean"],
                    "note": (
                        f"true {r['true_value']:.3g}; 95% HDI "
                        f"[{r['hdi_low']:.3g}, {r['hdi_high']:.3g}]; "
                        f"covered={bool(r['covered'])}"
                    ),
                }
            )
    return pd.DataFrame(rows), "ok"


def table_s3(store: ArtifactStore) -> Result:
    status = store.bayesian_status()
    if not status:
        return None, "no Bayesian run metadata"
    rows: List[Dict[str, Any]] = []
    for model in status.get("models_fitted", []):
        conv = store.convergence(model) or {}
        metrics = conv.get("metrics", {})
        rows.append(
            {
                "model": model,
                "engine": metrics.get("engine"),
                "chains": metrics.get("n_chains"),
                "draws_per_chain": metrics.get("n_draws_per_chain"),
                "parameters": metrics.get("n_parameters"),
                "max_rhat": metrics.get("max_rhat"),
                "min_ess_bulk": metrics.get("min_ess_bulk"),
                "min_ess_tail": metrics.get("min_ess_tail"),
                "divergences": metrics.get("divergences"),
                "passed": conv.get("passed"),
                "failures": "; ".join(conv.get("failures") or []),
            }
        )
    ppc = store.posterior_predictive() or {}
    if ppc.get("status") == "ok":
        rows.append(
            {
                "model": "posterior predictive check",
                "engine": "",
                "max_rhat": None,
                "failures": (
                    f"observed detection rate {ppc.get('observed_detection_rate'):.4f}; "
                    f"{ppc.get('n_replicates')} replicates"
                ),
            }
        )
    return pd.DataFrame(rows), "ok"


def table_s4(store: ArtifactStore) -> Result:
    selection = store.selection_results()
    lineages = store.lineage_results()
    status = store.evolution_status() or {}
    if (selection is None or not len(selection)) and (lineages is None or not len(lineages)):
        reason = status.get("reason") or "evolutionary analyses were not run"
        if status.get("missing_tools"):
            reason += " (missing tools: " + ", ".join(status["missing_tools"]) + ")"
        return None, reason
    frames = []
    if lineages is not None and len(lineages):
        lin = lineages.copy()
        lin["section"] = "alignment and phylogeny"
        frames.append(lin)
    if selection is not None and len(selection):
        sel = selection.copy()
        sel["section"] = "selection analysis"
        frames.append(sel)
    return pd.concat(frames, ignore_index=True), "ok"


def table_s5(store: ArtifactStore) -> Result:
    sensitivity = store.sensitivity()
    perm = store.connectivity_permutation()
    if (sensitivity is None or not len(sensitivity)) and (perm is None or not len(perm)):
        return None, "sensitivity analyses and null models are unavailable"
    frames = []
    if sensitivity is not None and len(sensitivity):
        s = sensitivity.copy()
        s["section"] = "habitat contrast sensitivity"
        frames.append(s)
    if perm is not None and len(perm):
        p = perm.copy()
        p["section"] = "connectivity permutation null"
        frames.append(p)
    return pd.concat(frames, ignore_index=True), "ok"


BUILDERS: Dict[str, Callable[[ArtifactStore], Result]] = {
    "Table1": table1,
    "Table2": table2,
    "Table3": table3,
    "Table4": table4,
    "S1": table_s1,
    "S2": table_s2,
    "S3": table_s3,
    "S4": table_s4,
    "S5": table_s5,
}


def build_tables(store: ArtifactStore, config: Dict[str, Any],
                 out_dir: str | Path) -> Dict[str, Any]:
    out = Path(out_dir)
    data_dir = out / "tables"
    md_dir = out / "manuscript"
    results: Dict[str, Any] = {}
    synthetic = store.is_synthetic()

    for name, spec in (config.get("tables") or {}).items():
        builder = BUILDERS.get(name)
        if builder is None:
            continue
        try:
            df, status = builder(store)
        except Exception as exc:
            df, status = None, f"generation failed: {type(exc).__name__}: {exc}"
        title = spec.get("title", name)
        if df is None or not len(df):
            results[name] = {"status": "not_generated", "reason": status, "title": title}
            write_text(
                f"# {name}. {title}\n\n**NOT GENERATED** — {status}\n",
                md_dir / f"{name}.md",
            )
            continue
        path = write_table(df, data_dir / name)
        banner = (
            "> **SYNTHETIC DATA — software validation only. Not a scientific result.**\n\n"
            if synthetic else ""
        )
        write_text(
            f"# {name}. {title}\n\n{banner}{to_markdown(df)}\n\n"
            f"_Rows: {len(df)}. Data file: `{Path(path).name}`._\n",
            md_dir / f"{name}.md",
        )
        results[name] = {"status": "ok", "title": title, "rows": int(len(df)),
                         "path": str(path)}
    return results
