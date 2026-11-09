"""Prior predictive simulation and parameter recovery."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from evonexus.bayesian.design import ModelDesign, design_matrix
from evonexus.bayesian.diagnostics import hdi
from evonexus.bayesian.predict import expit
from evonexus.core.io import read_json


def _draw_prior(rng: np.random.Generator, spec: dict[str, Any], size: int = 1) -> np.ndarray:
    kind = str(spec.get("dist", "normal")).lower()
    if kind == "normal":
        return rng.normal(float(spec.get("mu", 0.0)), float(spec.get("sigma", 1.0)), size)
    if kind in {"halfnormal", "half_normal"}:
        return np.abs(rng.normal(0.0, float(spec.get("sigma", 1.0)), size))
    if kind == "beta":
        return rng.beta(float(spec.get("alpha", 2.0)), float(spec.get("beta", 2.0)), size)
    if kind in {"exponential", "exp"}:
        return rng.exponential(1.0 / float(spec.get("lam", 1.0)), size)
    raise ValueError(f"unsupported prior distribution {kind!r}")


def prior_predictive_check(
    design: ModelDesign,
    priors: dict[str, Any],
    draws: int = 1000,
    seed: int = 20250101,
    plausible_band: tuple[float, float] = (0.0005, 0.35),
) -> dict[str, Any]:
    """Simulate detections under the priors alone and judge plausibility.

    Priors that imply, say, a 90% chance that every library carries every
    lineage are not weakly informative — they are wrong, and this check is how
    that is noticed before any data is used.
    """
    rng = np.random.default_rng(seed)
    X, _names, spans = design_matrix(design)
    _n, p = X.shape
    rates = np.empty(draws)
    per_draw_eta_sd = np.empty(draws)

    for d in range(draws):
        beta = np.zeros(p)
        for label, start, stop in spans:
            if label == "intercept":
                beta[start:stop] = _draw_prior(rng, priors.get("intercept", {}), stop - start)
            elif label.startswith("fixed:"):
                name = label.split(":", 1)[1]
                key = design.fixed_prior_keys[design.fixed_names.index(name)]
                beta[start:stop] = _draw_prior(rng, priors.get(key, {}), stop - start)
            else:
                block_name = label.split(":", 1)[1]
                block = design.block(block_name)
                scale = float(_draw_prior(rng, priors.get(block.prior_key, {}), 1)[0])
                if block.covariance is not None and not np.allclose(block.covariance, np.eye(block.size)):
                    chol = np.linalg.cholesky(block.covariance + 1e-8 * np.eye(block.size))
                    beta[start:stop] = scale * (chol @ rng.standard_normal(block.size))
                else:
                    beta[start:stop] = scale * rng.standard_normal(stop - start)
        eta = X @ beta
        per_draw_eta_sd[d] = float(np.std(eta))
        rates[d] = float(np.mean(expit(eta)))

    lo, hi = plausible_band
    inside = float(np.mean((rates >= lo) & (rates <= hi)))
    q = np.quantile(rates, [0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995])
    warnings: list[str] = []
    if inside < 0.5:
        warnings.append(
            f"only {inside:.0%} of prior predictive detection rates fall inside the "
            f"plausible band {lo}-{hi}; the priors imply implausible virome prevalence"
        )
    if q[3] > hi:
        warnings.append(f"median prior predictive detection rate {q[3]:.3f} exceeds {hi}")

    return {
        "status": "ok",
        "n_draws": draws,
        "detection_rate_quantiles": {
            "0.5%": float(q[0]), "2.5%": float(q[1]), "25%": float(q[2]),
            "50%": float(q[3]), "75%": float(q[4]), "97.5%": float(q[5]), "99.5%": float(q[6]),
        },
        "fraction_in_plausible_band": inside,
        "plausible_band": [lo, hi],
        "linear_predictor_sd_median": float(np.median(per_draw_eta_sd)),
        "warnings": warnings,
        "passed": not warnings,
        "samples": rates.tolist(),
    }


# ---------------------------------------------------------------------------
# parameter recovery
# ---------------------------------------------------------------------------
_RECOVERY_MAP = {
    "alpha": ("intercept", None),
    "beta_habitat": ("fixed", "habitat"),
    "gamma_depth": ("fixed", "depth"),
    "gamma_poolsize": ("fixed", "pool_size"),
    "gamma_time": ("fixed", "time"),
    "sigma_family": ("sigma", "virus_family"),
    "sigma_lineage": ("sigma", "virus_lineage"),
    "sigma_host": ("sigma", "host"),
    "sigma_tissue": ("sigma", "tissue"),
    "sigma_site": ("sigma", "site"),
    "sigma_interaction": ("sigma", "host_family_ix"),
}


def parameter_recovery(
    posterior: dict[str, Any],
    truth_path: str | Path,
    hdi_prob: float = 0.95,
) -> dict[str, Any]:
    """Compare posterior intervals with the parameters that generated the data.

    Only meaningful for synthetic data. The result is explicitly marked as
    synthetic so it can never be read as an empirical finding.
    """
    from evonexus.bayesian.predict import PosteriorPredictor

    path = Path(truth_path)
    if not path.exists():
        return {"status": "unavailable", "reason": f"no truth file at {path}"}
    truth_doc = read_json(path)
    truth = truth_doc.get("parameters", {})
    predictor = PosteriorPredictor(posterior)

    rows: list[dict[str, Any]] = []
    for name, value in truth.items():
        target = _RECOVERY_MAP.get(name)
        if target is None:
            continue
        kind, key = target
        if kind == "intercept":
            draws = predictor.intercept()
        elif kind == "fixed":
            if key not in predictor.fixed_span:
                continue
            draws = predictor.fixed(key)
        else:
            if key not in predictor.sigma_names:
                continue
            draws = predictor.sigma[:, predictor.sigma_names.index(key)]
        lo, hi = hdi(draws, hdi_prob)
        rows.append(
            {
                "parameter": name,
                "true_value": float(value),
                "posterior_mean": float(np.mean(draws)),
                "posterior_sd": float(np.std(draws, ddof=1)),
                "hdi_low": lo,
                "hdi_high": hi,
                "covered": bool(lo <= float(value) <= hi),
                "bias": float(np.mean(draws) - float(value)),
                "z_score": float((np.mean(draws) - float(value)) / max(np.std(draws, ddof=1), 1e-12)),
            }
        )

    table = pd.DataFrame(rows)
    coverage = float(table["covered"].mean()) if len(table) else float("nan")
    return {
        "status": "ok",
        "synthetic": True,
        "table": table,
        "coverage": coverage,
        "n_parameters_checked": len(table),
        "hdi_prob": hdi_prob,
        "note": (
            "Parameter recovery on SYNTHETIC data. It validates the software, "
            "not any biological claim."
        ),
    }
