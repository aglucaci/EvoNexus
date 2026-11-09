"""Convergence diagnostics and posterior summaries.

ArviZ is used when installed. The numpy implementations below are the
fallback, following the split-chain rank-normalized R-hat and effective sample
size definitions of Vehtari et al. (2021), so the numbers mean the same thing
whichever path runs.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats as sps


def _arviz():
    try:
        import arviz as az

        return az
    except Exception:
        return None


# ---------------------------------------------------------------------------
# numpy implementations
# ---------------------------------------------------------------------------
def _split_chains(x: np.ndarray) -> np.ndarray:
    """(chain, draw) -> (2*chain, draw/2)."""
    _n_chain, n_draw = x.shape
    half = n_draw // 2
    if half < 2:
        return x
    return np.concatenate([x[:, :half], x[:, half: 2 * half]], axis=0)


def _rank_normalize(x: np.ndarray) -> np.ndarray:
    flat = x.reshape(-1)
    ranks = sps.rankdata(flat, method="average")
    n = ranks.size
    z = sps.norm.ppf((ranks - 0.375) / (n + 0.25))
    return z.reshape(x.shape)


def rhat(x: np.ndarray) -> float:
    """Split rank-normalized R-hat for one parameter, shape (chain, draw)."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or x.shape[1] < 4:
        return float("nan")
    if np.allclose(x, x.flat[0]):
        return 1.0
    z = _rank_normalize(x)
    s = _split_chains(z)
    _m, n = s.shape
    chain_means = s.mean(axis=1)
    chain_vars = s.var(axis=1, ddof=1)
    W = chain_vars.mean()
    B = n * chain_means.var(ddof=1)
    if W <= 0:
        return float("nan")
    var_hat = ((n - 1) / n) * W + B / n
    return float(np.sqrt(var_hat / W))


def _autocov(series: np.ndarray) -> np.ndarray:
    n = series.size
    centred = series - series.mean()
    size = int(2 ** np.ceil(np.log2(2 * n)))
    freq = np.fft.rfft(centred, size)
    acov = np.fft.irfft(freq * np.conjugate(freq), size)[:n].real
    return acov / n


def ess(x: np.ndarray) -> float:
    """Bulk effective sample size for one parameter, shape (chain, draw)."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or x.shape[1] < 8:
        return float("nan")
    if np.allclose(x, x.flat[0]):
        return float(x.size)
    z = _rank_normalize(x)
    s = _split_chains(z)
    m, n = s.shape
    acov = np.stack([_autocov(chain) for chain in s])
    chain_mean = s.mean(axis=1)
    mean_var = acov[:, 0].mean() * n / (n - 1)
    var_plus = mean_var * (n - 1) / n
    if m > 1:
        var_plus += chain_mean.var(ddof=1)
    rho = np.zeros(n)
    rho[0] = 1.0
    if var_plus <= 0:
        return float("nan")
    for t in range(1, n):
        rho[t] = 1.0 - (mean_var - acov[:, t].mean()) / var_plus
    # Geyer initial positive sequence.
    tau = -1.0
    t = 0
    while t + 1 < n:
        pair = rho[t] + rho[t + 1]
        if pair < 0:
            break
        tau += 2 * pair
        t += 2
    tau = max(tau, 1.0 / np.log10(max(m * n, 11)))
    return float(m * n / tau)


def ess_tail(x: np.ndarray) -> float:
    """Tail ESS: the minimum of the 5% and 95% indicator ESS."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 2:
        return float("nan")
    flat = x.reshape(-1)
    lo, hi = np.quantile(flat, [0.05, 0.95])
    return float(min(ess((x <= lo).astype(float)), ess((x <= hi).astype(float))))


def mcse_mean(x: np.ndarray) -> float:
    n_eff = ess(x)
    if not np.isfinite(n_eff) or n_eff <= 0:
        return float("nan")
    return float(np.std(x, ddof=1) / np.sqrt(n_eff))


# ---------------------------------------------------------------------------
# summaries
# ---------------------------------------------------------------------------
def hdi(samples: np.ndarray, prob: float = 0.95) -> tuple[float, float]:
    """Highest-density interval of a flat sample vector."""
    s = np.sort(np.asarray(samples, dtype=float).reshape(-1))
    n = s.size
    if n == 0:
        return (float("nan"), float("nan"))
    k = max(1, int(np.floor(prob * n)))
    if k >= n:
        return (float(s[0]), float(s[-1]))
    widths = s[k:] - s[: n - k]
    i = int(np.argmin(widths))
    return (float(s[i]), float(s[i + k]))


def summarize_posterior(
    draws: np.ndarray,
    names: Sequence[str],
    hdi_prob: float = 0.95,
    max_parameters: int | None = None,
) -> pd.DataFrame:
    """Per-parameter summary table from an array of shape (chain, draw, p)."""
    arr = np.asarray(draws, dtype=float)
    if arr.ndim == 2:
        arr = arr[None, ...]
    n_params = arr.shape[-1]
    idx = range(n_params if max_parameters is None else min(n_params, max_parameters))
    rows: list[dict[str, Any]] = []
    for j in idx:
        x = arr[..., j]
        flat = x.reshape(-1)
        lo, hi = hdi(flat, hdi_prob)
        rows.append(
            {
                "parameter": names[j] if j < len(names) else f"param_{j}",
                "mean": float(np.mean(flat)),
                "sd": float(np.std(flat, ddof=1)) if flat.size > 1 else float("nan"),
                "hdi_low": lo,
                "hdi_high": hi,
                f"q{int((1 - hdi_prob) / 2 * 100)}": float(np.quantile(flat, (1 - hdi_prob) / 2)),
                "median": float(np.median(flat)),
                f"q{int((1 - (1 - hdi_prob) / 2) * 100)}": float(np.quantile(flat, 1 - (1 - hdi_prob) / 2)),
                "mcse_mean": mcse_mean(x),
                "ess_bulk": ess(x),
                "ess_tail": ess_tail(x),
                "r_hat": rhat(x),
            }
        )
    return pd.DataFrame(rows)


def convergence_report(posterior: dict[str, Any], gates: dict[str, Any]) -> dict[str, Any]:
    """Apply the convergence gates from ``config/priors.yaml``."""
    beta = np.asarray(posterior["beta"], dtype=float)
    names = list(posterior.get("coefficient_names") or [])
    summary = summarize_posterior(beta, names)

    sigma = posterior.get("sigma")
    if sigma is not None and np.asarray(sigma).size:
        sigma_summary = summarize_posterior(
            np.asarray(sigma, dtype=float),
            [f"sigma_{n}" for n in posterior.get("sigma_names", [])],
        )
        summary = pd.concat([summary, sigma_summary], ignore_index=True)

    stats = posterior.get("sampler_stats", {}) or {}
    n_total = int(beta.shape[0] * beta.shape[1])
    divergences = int(stats.get("divergences", 0))
    treedepth_hits = int(stats.get("treedepth_hits", 0))

    metrics = {
        "engine": posterior.get("engine"),
        "n_chains": int(beta.shape[0]),
        "n_draws_per_chain": int(beta.shape[1]),
        "n_parameters": int(beta.shape[2]),
        "max_rhat": float(np.nanmax(summary["r_hat"])) if len(summary) else float("nan"),
        "min_ess_bulk": float(np.nanmin(summary["ess_bulk"])) if len(summary) else float("nan"),
        "min_ess_tail": float(np.nanmin(summary["ess_tail"])) if len(summary) else float("nan"),
        "divergences": divergences,
        "divergence_fraction": divergences / max(1, n_total),
        "treedepth_hit_fraction": treedepth_hits / max(1, n_total),
    }

    failures: list[str] = []
    if np.isfinite(metrics["max_rhat"]) and metrics["max_rhat"] > float(gates.get("max_rhat", 1.01)):
        worst = summary.loc[summary["r_hat"].idxmax(), "parameter"]
        failures.append(
            f"max R-hat {metrics['max_rhat']:.3f} exceeds {gates.get('max_rhat')} "
            f"(worst parameter: {worst})"
        )
    if np.isfinite(metrics["min_ess_bulk"]) and metrics["min_ess_bulk"] < float(gates.get("min_ess_bulk", 400)):
        failures.append(
            f"minimum bulk ESS {metrics['min_ess_bulk']:.0f} below {gates.get('min_ess_bulk')}"
        )
    if np.isfinite(metrics["min_ess_tail"]) and metrics["min_ess_tail"] < float(gates.get("min_ess_tail", 400)):
        failures.append(
            f"minimum tail ESS {metrics['min_ess_tail']:.0f} below {gates.get('min_ess_tail')}"
        )
    if metrics["divergence_fraction"] > float(gates.get("max_divergences_fraction", 0.005)):
        failures.append(
            f"{divergences} divergent transitions "
            f"({metrics['divergence_fraction']:.2%} of draws)"
        )
    if metrics["treedepth_hit_fraction"] > float(gates.get("max_treedepth_hits_fraction", 0.05)):
        failures.append(
            f"{metrics['treedepth_hit_fraction']:.2%} of transitions hit the maximum tree depth"
        )

    return {
        "metrics": metrics,
        "failures": failures,
        "passed": not failures,
        "summary": summary,
        "note": stats.get("divergence_note"),
    }


# ---------------------------------------------------------------------------
# information criteria
# ---------------------------------------------------------------------------
def psis_loo(log_lik: np.ndarray) -> dict[str, Any]:
    """PSIS-LOO via ArviZ when available; otherwise plain importance-sampled LOO.

    The fallback reports ``method: "is_loo"`` and no Pareto-k diagnostics, so a
    comparison computed without ArviZ is never presented as if it had them.
    """
    az = _arviz()
    arr = np.asarray(log_lik, dtype=float)
    if arr.ndim == 2:
        arr = arr[None, ...]
    _n_chain, _n_draw, n_obs = arr.shape

    if az is not None:
        try:

            idata = az.from_dict(
                log_likelihood={"y": arr},
            )
            res = az.loo(idata, pointwise=True)
            return {
                "method": "psis_loo",
                "elpd_loo": float(res.elpd_loo),
                "p_loo": float(res.p_loo),
                "se": float(res.se),
                "pareto_k_gt_0_7": int(np.sum(np.asarray(res.pareto_k) > 0.7)),
                "n_obs": n_obs,
                "pointwise": np.asarray(res.loo_i).tolist() if hasattr(res, "loo_i") else None,
            }
        except Exception:
            pass

    flat = arr.reshape(-1, n_obs)
    s = flat.shape[0]
    # log of the harmonic-mean style importance weights for LOO.
    neg = -flat
    max_neg = neg.max(axis=0)
    lse = max_neg + np.log(np.exp(neg - max_neg).sum(axis=0))
    elpd_i = np.log(s) - lse
    lppd_i = -np.log(s) + (
        flat.max(axis=0) + np.log(np.exp(flat - flat.max(axis=0)).sum(axis=0))
    )
    return {
        "method": "is_loo",
        "elpd_loo": float(elpd_i.sum()),
        "p_loo": float((lppd_i - elpd_i).sum()),
        "se": float(np.sqrt(n_obs) * np.std(elpd_i, ddof=1)),
        "pareto_k_gt_0_7": None,
        "n_obs": n_obs,
        "pointwise": elpd_i.tolist(),
        "caveat": (
            "ArviZ is not installed; this is plain importance-sampled LOO without "
            "Pareto smoothing or k diagnostics. Install ArviZ for published results."
        ),
    }


def compare_models(loos: dict[str, dict[str, Any]]) -> pd.DataFrame:
    """Rank models by ELPD with differences and their standard errors."""
    rows = []
    best_name = max(loos, key=lambda k: loos[k].get("elpd_loo", -np.inf))
    best = np.asarray(loos[best_name].get("pointwise") or [], dtype=float)
    for name, loo in loos.items():
        pointwise = np.asarray(loo.get("pointwise") or [], dtype=float)
        if pointwise.size and best.size == pointwise.size:
            diff = pointwise - best
            elpd_diff = float(diff.sum())
            dse = float(np.sqrt(diff.size) * np.std(diff, ddof=1)) if diff.size > 1 else float("nan")
        else:
            elpd_diff = float(loo.get("elpd_loo", np.nan) - loos[best_name].get("elpd_loo", np.nan))
            dse = float("nan")
        rows.append(
            {
                "model": name,
                "method": loo.get("method"),
                "elpd_loo": loo.get("elpd_loo"),
                "p_loo": loo.get("p_loo"),
                "se": loo.get("se"),
                "elpd_diff": elpd_diff,
                "dse": dse,
                "pareto_k_gt_0_7": loo.get("pareto_k_gt_0_7"),
            }
        )
    out = pd.DataFrame(rows).sort_values("elpd_loo", ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out
