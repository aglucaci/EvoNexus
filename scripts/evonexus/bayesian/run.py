"""Driver for the Bayesian inference stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from evonexus.bayesian.checks import parameter_recovery, prior_predictive_check
from evonexus.bayesian.contrast import (
    ESTIMAND_CAVEAT,
    habitat_contrast,
    posterior_predictive_replicates,
    raw_detection_rates,
)
from evonexus.bayesian.design import build_design
from evonexus.bayesian.diagnostics import compare_models, convergence_report, psis_loo
from evonexus.bayesian.fit import fit_model, resolve_engine, save_posterior
from evonexus.core.config import load_config, resolve_variant
from evonexus.core.io import read_json, read_table, write_json, write_table
from evonexus.core.provenance import RunRecord
from evonexus.quantify.matrix import apply_restriction


def _host_covariance(path: str | Path | None) -> Dict[str, Any] | None:
    """Load a host phylogenetic covariance matrix if one was produced."""
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return None
    doc = read_json(p)
    if "species" in doc and "matrix" in doc:
        return doc
    return None


def _model_terms(comparisons: Dict[str, Any], model: str) -> List[str]:
    spec = (comparisons.get("models") or {}).get(model)
    if spec is None:
        raise KeyError(f"model {model!r} is not defined in model_comparisons.yaml")
    return list(spec.get("terms") or [])


def _requirements_met(spec: Dict[str, Any], context: Dict[str, bool]) -> tuple[bool, List[str]]:
    unmet = [k for k, v in (spec.get("requires") or {}).items() if bool(v) and not context.get(k, False)]
    return (not unmet), unmet


def run_bayesian(
    frame_path: str | Path,
    out_dir: str | Path,
    models: List[str] | None = None,
    priors_config: str = "priors",
    comparisons_config: str = "model_comparisons",
    prior_variant: str = "regularizing",
    engine: str = "auto",
    draws: int | None = None,
    tune: int | None = None,
    chains: int | None = None,
    seed: int | None = None,
    host_covariance_path: str | Path | None = None,
    truth_path: str | Path | None = None,
    synthetic: bool = False,
    run_sensitivity: bool = True,
    skip_primary: bool = False,
    processed_dir: str | Path | None = None,
) -> Dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    priors_cfg = load_config(priors_config)
    comparisons = load_config(comparisons_config)
    priors = resolve_variant(priors_cfg, prior_variant)
    sampler = dict(priors_cfg.get("sampler") or {})
    for key, value in (("draws", draws), ("tune", tune), ("chains", chains), ("seed", seed)):
        if value is not None:
            sampler[key] = value
    gates = priors_cfg.get("convergence_gates") or {}

    record = RunRecord(
        step="bayesian",
        params={
            "models": models, "prior_variant": prior_variant, "engine": engine,
            "sampler": sampler, "synthetic": synthetic, "skip_primary": skip_primary,
        },
    )
    record.synthetic = synthetic
    record.add_input(frame_path, role="analysis_frame", checksum=False)

    frame = read_table(frame_path)
    if frame.empty:
        record.error("analysis frame is empty; no model was fitted")
        record.finish(out / "run_record_bayesian.json", status="empty")
        return {"status": "empty", "reason": "analysis frame is empty"}

    if skip_primary:
        # A critical QC gate failed. Descriptive output only — and it says so.
        descriptive = raw_detection_rates(frame)
        path = write_table(descriptive, out / "descriptive_detection_rates")
        write_json(
            {
                "status": "skipped",
                "reason": "a critical QC gate failed; the primary habitat model was not run",
                "descriptive_only": True,
                "message": (
                    "No habitat effect is estimated or implied. The table lists "
                    "unadjusted per-library detection rates only."
                ),
            },
            out / "bayesian_status.json",
        )
        record.warn("primary model skipped because of a critical QC gate failure")
        record.add_output(path, role="descriptive")
        record.finish(out / "run_record_bayesian.json", status="skipped")
        return {"status": "skipped", "descriptive": str(path)}

    chosen_engine = resolve_engine(engine)
    record.note(f"sampler engine: {chosen_engine}")
    if chosen_engine != "pymc_nuts":
        record.warn(
            "PyMC is not installed; using the built-in Polya-Gamma Gibbs engine. "
            "Published results should be produced with PyMC/NUTS."
        )

    ecology_columns = [c for c in frame.columns if c.endswith("_z")
                       and c not in {"log_depth_z", "log_pool_size_z", "time_z"}]
    host_cov = _host_covariance(host_covariance_path)
    context = {
        "ecology_covariates_available": bool(ecology_columns),
        "host_phylogeny_available": host_cov is not None,
    }

    models = models or list((comparisons.get("models") or {}).keys())
    primary = comparisons.get("primary_model", "M3")

    results: Dict[str, Any] = {}
    loos: Dict[str, Dict[str, Any]] = {}
    fitted: Dict[str, Dict[str, Any]] = {}
    skipped: Dict[str, List[str]] = {}

    prior_check_written = False

    for model_name in models:
        spec = (comparisons.get("models") or {}).get(model_name) or {}
        ok, unmet = _requirements_met(spec, context)
        if not ok:
            skipped[model_name] = unmet
            record.warn(f"{model_name} not fitted: unmet requirements {unmet}")
            continue

        terms = _model_terms(comparisons, model_name)
        design = build_design(
            frame, terms, model_name=model_name,
            ecology_columns=ecology_columns, host_covariance=host_cov,
        )
        record.note(
            f"{model_name}: {design.n_obs} observations, "
            f"{design.parameter_count()} parameters, blocks={design.summary()['blocks']}"
        )

        if not prior_check_written:
            band = ((priors_cfg.get("prior_predictive") or {}).get("plausible_detection_rate")
                    or [0.0005, 0.35])
            pc = prior_predictive_check(
                design, priors,
                draws=int((priors_cfg.get("prior_predictive") or {}).get("draws", 1000)),
                seed=int(sampler.get("seed", 20250101)),
                plausible_band=(float(band[0]), float(band[1])),
            )
            samples = pc.pop("samples")
            pc["n_samples_recorded"] = len(samples)
            write_json(pc, out / "prior_predictive_check.json")
            write_table(pd.DataFrame({"prior_predictive_detection_rate": samples}),
                        out / "prior_predictive_draws")
            for w in pc["warnings"]:
                record.warn(f"prior predictive: {w}")
            prior_check_written = True

        posterior = fit_model(design, priors, sampler=sampler, engine=engine)
        paths = save_posterior(posterior, design, out, model_name)
        for p in paths.values():
            record.add_output(p, role=f"posterior:{model_name}", checksum=False)

        conv = convergence_report(posterior, gates)
        write_table(conv["summary"], out / f"posterior_summary_{model_name}")
        write_json(
            {k: v for k, v in conv.items() if k != "summary"},
            out / f"convergence_{model_name}.json",
        )
        if not conv["passed"]:
            for f in conv["failures"]:
                record.warn(f"{model_name} convergence: {f}")

        if posterior.get("log_likelihood") is not None:
            loos[model_name] = psis_loo(posterior["log_likelihood"])

        fitted[model_name] = posterior
        results[model_name] = {
            "paths": paths,
            "convergence": conv["metrics"],
            "convergence_passed": conv["passed"],
            "convergence_failures": conv["failures"],
            "design": design.summary(),
            "design_metadata": design.metadata,
        }

    if not fitted:
        write_json({"status": "no_models_fitted", "skipped": skipped}, out / "bayesian_status.json")
        record.error("no model could be fitted")
        record.finish(out / "run_record_bayesian.json", status="failed")
        return {"status": "no_models_fitted", "skipped": skipped}

    comparison_path = None
    if len(loos) > 1:
        comparison = compare_models(loos)
        comparison_path = write_table(comparison, out / "model_comparison")
        write_json({k: {kk: vv for kk, vv in v.items() if kk != "pointwise"}
                    for k, v in loos.items()}, out / "model_loo.json")
        record.add_output(comparison_path, role="model_comparison")

    # ---- primary estimand ------------------------------------------------
    primary_name = primary if primary in fitted else sorted(fitted)[-1]
    primary_post = fitted[primary_name]
    contrast = habitat_contrast(primary_post, frame, seed=int(sampler.get("seed", 20250101)))

    contrast_paths: Dict[str, str] = {}
    if contrast.get("status") == "ok":
        contrast_paths["table"] = str(write_table(contrast["table"], out / "habitat_contrast"))
        contrast_paths["by_family"] = str(write_table(contrast["by_family"], out / "habitat_contrast_by_family"))
        contrast_paths["by_host"] = str(write_table(contrast["by_host"], out / "habitat_contrast_by_host"))
        write_table(
            pd.DataFrame(
                {
                    "p_settlement": contrast["draws"]["p_settlement"],
                    "p_wild": contrast["draws"]["p_wild"],
                    "difference": contrast["draws"]["difference"],
                }
            ),
            out / "habitat_contrast_draws",
        )
        write_json(
            {
                "status": "ok",
                "model": primary_name,
                "engine": contrast["engine"],
                "reference": contrast["reference"],
                "caveat": ESTIMAND_CAVEAT,
                "synthetic": synthetic,
            },
            out / "habitat_contrast_metadata.json",
        )
        row = contrast["table"].iloc[2]
        record.note(
            f"standardized difference (settlement - wild): "
            f"{row['mean']:.4f} [{row['hdi_95_low']:.4f}, {row['hdi_95_high']:.4f}] — "
            f"{ESTIMAND_CAVEAT}"
        )
    else:
        record.warn(f"habitat contrast not computed: {contrast.get('reason')}")

    write_table(raw_detection_rates(frame), out / "raw_detection_rates")

    ppc = posterior_predictive_replicates(primary_post, frame,
                                          seed=int(sampler.get("seed", 20250101)))
    write_json(ppc, out / "posterior_predictive_check.json")

    recovery = None
    if truth_path:
        recovery = parameter_recovery(primary_post, truth_path)
        if recovery.get("status") == "ok":
            write_table(recovery["table"], out / "parameter_recovery")
            write_json({k: v for k, v in recovery.items() if k != "table"},
                       out / "parameter_recovery.json")
            record.note(
                f"parameter recovery coverage: {recovery['coverage']:.2f} "
                f"over {recovery['n_parameters_checked']} parameters (synthetic)"
            )

    # ---- sensitivity -----------------------------------------------------
    sensitivity_path = None
    if run_sensitivity and contrast.get("status") == "ok":
        sens_rows: List[Dict[str, Any]] = []
        sens_cfg = comparisons.get("sensitivity") or {}
        tables = None
        if processed_dir:
            from evonexus.qc.run import load_canonical

            tables = load_canonical(processed_dir)

        for variant in sens_cfg.get("prior_variants", []):
            if variant == prior_variant:
                sens_rows.append(_sens_row("prior", variant, contrast))
                continue
            try:
                alt_priors = resolve_variant(priors_cfg, variant)
                design = build_design(frame, _model_terms(comparisons, primary_name),
                                      model_name=f"{primary_name}_{variant}",
                                      ecology_columns=ecology_columns, host_covariance=host_cov)
                alt_sampler = dict(sampler)
                alt_sampler["draws"] = max(200, int(sampler.get("draws", 1000)) // 2)
                alt_post = fit_model(design, alt_priors, sampler=alt_sampler, engine=engine,
                                     store_log_lik=False)
                alt_contrast = habitat_contrast(alt_post, frame,
                                                seed=int(sampler.get("seed", 20250101)))
                sens_rows.append(_sens_row("prior", variant, alt_contrast))
            except Exception as exc:
                record.warn(f"sensitivity prior={variant} failed: {type(exc).__name__}: {exc}")

        for restriction in (sens_cfg.get("restrictions") or []):
            name = restriction.get("name")
            try:
                sub = apply_restriction(frame, name, tables)
                if len(sub) < 50 or sub["habitat"].nunique() < 2:
                    sens_rows.append(
                        {"dimension": "restriction", "variant": name, "status": "not_estimable",
                         "reason": f"{len(sub)} cells / {sub['habitat'].nunique()} habitats remain"}
                    )
                    continue
                design = build_design(sub, _model_terms(comparisons, primary_name),
                                      model_name=f"{primary_name}_{name}",
                                      ecology_columns=[c for c in ecology_columns if c in sub.columns],
                                      host_covariance=host_cov)
                alt_sampler = dict(sampler)
                alt_sampler["draws"] = max(200, int(sampler.get("draws", 1000)) // 2)
                alt_post = fit_model(design, priors, sampler=alt_sampler, engine=engine,
                                     store_log_lik=False)
                alt_contrast = habitat_contrast(alt_post, sub,
                                                seed=int(sampler.get("seed", 20250101)))
                sens_rows.append(_sens_row("restriction", name, alt_contrast, n_cells=len(sub)))
            except Exception as exc:
                record.warn(f"sensitivity restriction={name} failed: {type(exc).__name__}: {exc}")

        if sens_rows:
            sensitivity_path = write_table(pd.DataFrame(sens_rows), out / "sensitivity")
            record.add_output(sensitivity_path, role="sensitivity")

    status_payload = {
        "status": "ok",
        "engine": chosen_engine,
        "primary_model": primary_name,
        "models_fitted": sorted(fitted),
        "models_skipped": skipped,
        "prior_variant": prior_variant,
        "sampler": sampler,
        "synthetic": synthetic,
        "caveat": ESTIMAND_CAVEAT,
        "convergence": {k: v["convergence_passed"] for k, v in results.items()},
    }
    write_json(status_payload, out / "bayesian_status.json")
    all_converged = all(v["convergence_passed"] for v in results.values())
    record.finish(out / "run_record_bayesian.json",
                  status="ok" if all_converged else "convergence_warning")

    return {
        "status": "ok",
        "primary_model": primary_name,
        "models": results,
        "comparison": str(comparison_path) if comparison_path else None,
        "contrast": contrast_paths,
        "sensitivity": str(sensitivity_path) if sensitivity_path else None,
        "recovery": recovery["status"] if recovery else None,
        "engine": chosen_engine,
    }


def _sens_row(dimension: str, variant: str, contrast: Dict[str, Any],
              n_cells: int | None = None) -> Dict[str, Any]:
    if contrast.get("status") != "ok":
        return {"dimension": dimension, "variant": variant, "status": "failed",
                "reason": contrast.get("reason")}
    row = contrast["table"].iloc[2]
    return {
        "dimension": dimension,
        "variant": variant,
        "status": "ok",
        "n_cells": n_cells,
        "mean_difference": float(row["mean"]),
        "hdi_low": float(row["hdi_95_low"]),
        "hdi_high": float(row["hdi_95_high"]),
        "p_direction_positive": float(row["p_direction_positive"]),
    }
