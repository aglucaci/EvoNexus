"""Held-out-site cross-validation and calibration.

The rule that makes this honest: when a site (or region) is held out, its
random effect is **integrated over its estimated population distribution**,
not estimated from the held-out observations. Estimating it from the data
being predicted would leak the answer and overstate generalization.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from evonexus.bayesian.design import build_design
from evonexus.bayesian.fit import fit_model
from evonexus.bayesian.predict import PosteriorPredictor, expit
from evonexus.core.config import load_config, resolve_variant
from evonexus.core.io import read_table, write_json, write_table
from evonexus.core.provenance import RunRecord


def scoring_rules(y: np.ndarray, p_draws: np.ndarray) -> dict[str, float]:
    """Log predictive score, Brier score and AUC from posterior draws of p."""
    y = np.asarray(y, dtype=float)
    p_mean = p_draws.mean(axis=0)
    eps = 1e-12
    # Pointwise log predictive density, averaging over the posterior first.
    lpd = np.log(np.clip(np.where(y > 0.5, p_draws, 1.0 - p_draws).mean(axis=0), eps, 1.0))
    brier = float(np.mean((p_mean - y) ** 2))

    order = np.argsort(p_mean)
    n_pos, n_neg = float(y.sum()), float((1 - y).sum())
    if n_pos > 0 and n_neg > 0:
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, len(y) + 1)
        auc = float((ranks[y > 0.5].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
    else:
        auc = float("nan")

    return {
        "elpd": float(lpd.sum()),
        "mean_lpd": float(lpd.mean()),
        "brier": brier,
        "auc": auc,
        "n": len(y),
        "observed_rate": float(y.mean()),
        "predicted_rate": float(p_mean.mean()),
    }


def calibration_curve(y: np.ndarray, p: np.ndarray, n_bins: int = 10) -> pd.DataFrame:
    """Observed frequency against predicted probability, in equal-count bins."""
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    if len(y) == 0:
        return pd.DataFrame()
    quantiles = np.quantile(p, np.linspace(0, 1, n_bins + 1))
    quantiles[0] -= 1e-9
    quantiles[-1] += 1e-9
    bins = np.unique(quantiles)
    idx = np.clip(np.digitize(p, bins) - 1, 0, len(bins) - 2)
    rows: list[dict[str, Any]] = []
    for b in range(len(bins) - 1):
        mask = idx == b
        if not mask.any():
            continue
        n = int(mask.sum())
        observed = float(y[mask].mean())
        se = float(np.sqrt(max(observed * (1 - observed), 1e-12) / n))
        rows.append(
            {
                "bin": b,
                "predicted_mean": float(p[mask].mean()),
                "observed_frequency": observed,
                "n": n,
                "se": se,
                "lower": max(0.0, observed - 1.96 * se),
                "upper": min(1.0, observed + 1.96 * se),
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        out["calibration_error"] = out["observed_frequency"] - out["predicted_mean"]
    return out


def _folds(frame: pd.DataFrame, group_column: str, n_folds: int,
           seed: int) -> list[tuple[str, np.ndarray]]:
    groups = sorted(set(frame[group_column].astype(str)))
    rng = np.random.default_rng(seed)
    rng.shuffle(groups)
    n_folds = min(n_folds, len(groups))
    if n_folds < 2:
        return []
    chunks = np.array_split(np.asarray(groups, dtype=object), n_folds)
    return [(f"fold{i+1}", np.asarray(chunk, dtype=object)) for i, chunk in enumerate(chunks)]


def cross_validate(
    frame: pd.DataFrame,
    terms: Sequence[str],
    priors: dict[str, Any],
    sampler: dict[str, Any],
    group_column: str = "site_id",
    n_folds: int = 5,
    engine: str = "auto",
    seed: int = 20250101,
    ecology_columns: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Refit on training groups, predict held-out groups with integration."""
    folds = _folds(frame, group_column, n_folds, seed)
    if not folds:
        return {"status": "not_estimable",
                "reason": f"fewer than two distinct {group_column} values"}

    fold_rows: list[dict[str, Any]] = []
    pooled_y: list[np.ndarray] = []
    pooled_p: list[np.ndarray] = []
    failures: list[str] = []

    for fold_name, held_out in folds:
        test_mask = frame[group_column].astype(str).isin(set(map(str, held_out)))
        train = frame[~test_mask]
        test = frame[test_mask]
        if train.empty or test.empty or train["habitat"].nunique() < 2:
            failures.append(f"{fold_name}: not estimable (train or test empty, or one habitat)")
            continue
        try:
            design = build_design(train, list(terms), model_name=f"cv_{fold_name}",
                                 ecology_columns=ecology_columns)
            posterior = fit_model(design, priors, sampler=sampler, engine=engine,
                                  store_log_lik=False)
            predictor = PosteriorPredictor(posterior)
            # Held-out site AND batch effects are integrated over, never fitted.
            p = expit(
                predictor.eta(test, integrate=(group_column.replace("site_id", "site"), "batch"),
                              rng=np.random.default_rng(seed))
            )
            y = pd.to_numeric(test["y"], errors="coerce").to_numpy(float)
            scores = scoring_rules(y, p)
            scores.update({"fold": fold_name, "held_out_groups": list(map(str, held_out)),
                           "n_train": len(train), "group_column": group_column})
            fold_rows.append(scores)
            pooled_y.append(y)
            pooled_p.append(p.mean(axis=0))
        except Exception as exc:
            failures.append(f"{fold_name}: {type(exc).__name__}: {exc}")

    if not fold_rows:
        return {"status": "failed", "failures": failures}

    y_all = np.concatenate(pooled_y)
    p_all = np.concatenate(pooled_p)
    return {
        "status": "ok",
        "folds": pd.DataFrame(fold_rows),
        "calibration": calibration_curve(y_all, p_all),
        "overall": {
            "elpd": float(sum(r["elpd"] for r in fold_rows)),
            "mean_lpd": float(np.mean([r["mean_lpd"] for r in fold_rows])),
            "brier": float(np.mean([r["brier"] for r in fold_rows])),
            "auc": float(np.nanmean([r["auc"] for r in fold_rows])),
            "n": len(y_all),
        },
        "failures": failures,
        "note": (
            "Held-out group random effects were integrated over their estimated "
            "population distribution, not estimated from held-out observations."
        ),
    }


def run_validation(
    frame_path: str | Path,
    out_dir: str | Path,
    posterior_dir: str | Path | None = None,
    model: str = "M3",
    group_columns: Sequence[str] = ("site_id",),
    n_folds: int = 5,
    priors_config: str = "priors",
    comparisons_config: str = "model_comparisons",
    prior_variant: str = "regularizing",
    engine: str = "auto",
    draws: int | None = None,
    tune: int | None = None,
    chains: int | None = None,
    synthetic: bool = False,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    record = RunRecord(step="validation",
                       params={"model": model, "group_columns": list(group_columns),
                               "n_folds": n_folds, "synthetic": synthetic})
    record.synthetic = synthetic

    frame = read_table(frame_path)
    if frame.empty:
        record.error("analysis frame is empty")
        record.finish(out / "run_record_validation.json", status="empty")
        return {"status": "empty"}

    priors_cfg = load_config(priors_config)
    comparisons = load_config(comparisons_config)
    priors = resolve_variant(priors_cfg, prior_variant)
    sampler = dict(priors_cfg.get("sampler") or {})
    for key, value in (("draws", draws), ("tune", tune), ("chains", chains)):
        if value is not None:
            sampler[key] = value
    terms = ((comparisons.get("models") or {}).get(model) or {}).get("terms") or []
    ecology_columns = [c for c in frame.columns if c.endswith("_z")
                       and c not in {"log_depth_z", "log_pool_size_z", "time_z"}]

    results: dict[str, Any] = {}
    all_folds: list[pd.DataFrame] = []
    for group in group_columns:
        if group not in frame.columns:
            record.warn(f"grouping column {group!r} absent; skipped")
            continue
        res = cross_validate(frame, terms, priors, sampler, group_column=group,
                             n_folds=n_folds, engine=engine,
                             ecology_columns=ecology_columns)
        results[group] = {k: v for k, v in res.items() if k not in {"folds", "calibration"}}
        if res.get("status") == "ok":
            folds = res["folds"].copy()
            folds["grouping"] = group
            all_folds.append(folds)
            write_table(res["calibration"], out / f"calibration_{group}")
            record.note(
                f"held-out {group}: mean lpd {res['overall']['mean_lpd']:.4f}, "
                f"Brier {res['overall']['brier']:.4f}, AUC {res['overall']['auc']:.3f}"
            )
        else:
            record.warn(f"cross-validation on {group} not completed: {res.get('status')}")

    folds_path = None
    if all_folds:
        folds_path = write_table(pd.concat(all_folds, ignore_index=True), out / "cv_folds")
        record.add_output(folds_path, role="cross_validation", checksum=False)

    # In-sample calibration of the fitted primary model, for comparison.
    insample = None
    if posterior_dir:
        from evonexus.bayesian.fit import load_posterior

        try:
            posterior = load_posterior(posterior_dir, model)
            predictor = PosteriorPredictor(posterior)
            p = expit(predictor.eta(frame, integrate=(), rng=np.random.default_rng(1)))
            y = pd.to_numeric(frame["y"], errors="coerce").to_numpy(float)
            insample = scoring_rules(y, p)
            write_table(calibration_curve(y, p.mean(axis=0)), out / "calibration_in_sample")
        except Exception as exc:
            record.warn(f"in-sample calibration unavailable: {type(exc).__name__}: {exc}")

    write_json(
        {
            "status": "ok" if all_folds else "not_estimable",
            "model": model,
            "cross_validation": results,
            "in_sample": insample,
            "synthetic": synthetic,
            "note": (
                "Cross-validated scores answer whether the model predicts detections "
                "at unseen locations. They do not validate the causal interpretation "
                "of any habitat difference."
            ),
        },
        out / "validation_metadata.json",
    )
    record.finish(out / "run_record_validation.json",
                  status="ok" if all_folds else "not_estimable")
    return {"status": "ok" if all_folds else "not_estimable",
            "folds": str(folds_path) if folds_path else None}
