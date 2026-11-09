"""A self-contained Gibbs sampler for the hierarchical detection model.

Why this exists
---------------
PyMC/NUTS is the intended engine for published results. But a reproducibility
pipeline that cannot run at all without a heavy dependency is not much use on
a laptop or in CI, so EvoNexus ships a second engine for the *same*
:class:`~evonexus.bayesian.design.ModelDesign`.

Method
------
Bayesian logistic regression with Polya-Gamma data augmentation
(Polson, Scott & Windle 2013). Conditional on latent weights
:math:`\\omega_i \\sim \\mathrm{PG}(1, x_i^\\top\\beta)` the coefficient vector
has a closed-form Gaussian conditional, which makes a blocked Gibbs sweep
exact for the coefficients. Hierarchical scales get a Metropolis step on
:math:`\\log\\sigma` with the configured half-normal prior.

Honest limitations, recorded in the output metadata
---------------------------------------------------
* The Polya-Gamma draws use the truncated infinite-sum representation with a
  deterministic tail correction; this is an approximation, not an exact
  Devroye sampler.
* Mixing is slower per iteration than NUTS, so convergence must be checked,
  not assumed. The same R-hat/ESS gates apply.
* ``engine`` is recorded as ``polya_gamma_gibbs`` in every artifact so a
  reader can never mistake a fallback fit for a NUTS fit.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from evonexus.bayesian.design import ModelDesign, design_matrix

ENGINE = "polya_gamma_gibbs"
_PG_TRUNCATION = 200


def _pg_sample(rng: np.random.Generator, c: np.ndarray, truncation: int = _PG_TRUNCATION) -> np.ndarray:
    """Draw omega ~ PG(1, c) elementwise via the truncated sum representation."""
    c = np.abs(np.asarray(c, dtype=float))
    n = c.shape[0]
    k = np.arange(1, truncation + 1, dtype=float)
    denom = (k - 0.5) ** 2 + (c[:, None] ** 2) / (4.0 * math.pi ** 2)
    g = rng.exponential(size=(n, truncation))
    omega = (g / denom).sum(axis=1) / (2.0 * math.pi ** 2)
    # Deterministic mean correction for the discarded tail.
    tail = 1.0 / (2.0 * math.pi ** 2) * (1.0 / (k[-1] + 0.5) + 0.0)
    return omega + tail


def _half_normal_logpdf(x: np.ndarray | float, sigma: float) -> float:
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        return -np.inf
    return float(np.sum(-0.5 * (x / sigma) ** 2 - math.log(sigma) + 0.5 * math.log(2.0 / math.pi)))


def _fixed_prior_sd(prior: dict[str, Any]) -> float:
    return float(prior.get("sigma", 1.0))


def _fixed_prior_mu(prior: dict[str, Any]) -> float:
    return float(prior.get("mu", 0.0))


def _block_prior_scale(prior: dict[str, Any]) -> float:
    return float(prior.get("sigma", 1.0))


def _slice_sample_log(
    rng: np.random.Generator,
    log_density: Any,
    current: float,
    width: float = 1.0,
    max_steps: int = 12,
    max_shrink: int = 60,
) -> float:
    """Stepping-out slice sampler for a positive scale, on the log scale.

    A variance component whose posterior mass sits near zero makes a random-walk
    Metropolis step stick, because log(sigma) runs off towards minus infinity.
    Slice sampling has no step size to tune and handles that geometry, which
    matters here: several blocks (sequencing batch, tissue) legitimately have
    near-zero scales.
    """
    def ld(t: float) -> float:
        # log density in t = log(sigma), including the Jacobian.
        return log_density(math.exp(t)) + t

    x0 = math.log(max(current, 1e-12))
    f0 = ld(x0)
    if not np.isfinite(f0):
        return current
    y = f0 + math.log(rng.random() + 1e-300)

    left = x0 - width * rng.random()
    right = left + width
    j = int(max_steps * rng.random())
    k = max_steps - 1 - j
    while j > 0 and y < ld(left):
        left -= width
        j -= 1
    while k > 0 and y < ld(right):
        right += width
        k -= 1

    for _ in range(max_shrink):
        x1 = left + rng.random() * (right - left)
        if y < ld(x1):
            return float(math.exp(x1))
        if x1 < x0:
            left = x1
        else:
            right = x1
    return current


def _chol_solve(precision: np.ndarray, rhs: np.ndarray, rng: np.random.Generator,
                jitter: float = 1e-8) -> np.ndarray:
    """Draw from N(precision^-1 rhs, precision^-1) using a Cholesky factor."""
    p = precision.shape[0]
    for attempt in range(6):
        try:
            L = np.linalg.cholesky(precision + (jitter * 10 ** attempt) * np.eye(p))
            break
        except np.linalg.LinAlgError:
            continue
    else:  # pragma: no cover - extremely ill-conditioned
        L = np.linalg.cholesky(precision + 1e-2 * np.eye(p))
    mean = np.linalg.solve(L.T, np.linalg.solve(L, rhs))
    z = rng.standard_normal(p)
    return mean + np.linalg.solve(L.T, z)


def sample(
    design: ModelDesign,
    priors: dict[str, Any],
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    seed: int = 20250101,
    thin: int = 1,
    store_log_lik: bool = True,
    progress_every: int = 0,
) -> dict[str, Any]:
    """Run the Gibbs sampler and return posterior draws in a plain dict."""
    X, names, spans = design_matrix(design)
    y = design.y
    n, p = X.shape
    kappa = y - 0.5

    block_spans: list[tuple[str, int, int, Any]] = []
    fixed_mu = np.zeros(p)
    fixed_prec = np.zeros(p)

    for label, start, stop in spans:
        if label == "intercept":
            prior = priors.get("intercept", {})
            fixed_mu[start:stop] = _fixed_prior_mu(prior)
            fixed_prec[start:stop] = 1.0 / _fixed_prior_sd(prior) ** 2
        elif label.startswith("fixed:"):
            name = label.split(":", 1)[1]
            idx = design.fixed_names.index(name)
            prior = priors.get(design.fixed_prior_keys[idx], {})
            fixed_mu[start:stop] = _fixed_prior_mu(prior)
            fixed_prec[start:stop] = 1.0 / _fixed_prior_sd(prior) ** 2
        else:
            block_name = label.split(":", 1)[1]
            block = design.block(block_name)
            block_spans.append((block_name, start, stop, block))

    n_blocks = len(block_spans)
    block_prior_scales = np.array(
        [_block_prior_scale(priors.get(b[3].prior_key, {})) for b in block_spans]
    ) if n_blocks else np.zeros(0)

    # Precompute inverse covariances for phylogenetic blocks.
    block_inv_cov: dict[str, np.ndarray | None] = {}
    for name, start, stop, block in block_spans:
        if block.covariance is not None and not np.allclose(block.covariance, np.eye(block.size)):
            block_inv_cov[name] = np.linalg.inv(
                block.covariance + 1e-8 * np.eye(block.size)
            )
        else:
            block_inv_cov[name] = None

    chain_beta: list[np.ndarray] = []
    chain_sigma: list[np.ndarray] = []
    chain_loglik: list[np.ndarray] = []

    kept = max(1, draws // max(1, thin))

    for chain in range(chains):
        rng = np.random.default_rng(seed + 1000 * chain)
        beta = rng.normal(0.0, 0.1, p)
        beta[:1] = fixed_mu[:1]
        sigma = np.full(n_blocks, 0.5) if n_blocks else np.zeros(0)
        accepted = np.zeros(n_blocks)
        proposals = np.zeros(n_blocks)

        beta_draws = np.empty((kept, p))
        sigma_draws = np.empty((kept, max(n_blocks, 0)))
        loglik_draws = np.empty((kept, n)) if store_log_lik else None

        store = 0
        total = tune + draws
        for it in range(total):
            eta = X @ beta
            omega = _pg_sample(rng, eta)

            precision = (X * omega[:, None]).T @ X
            prior_precision = np.zeros((p, p))
            np.fill_diagonal(prior_precision, fixed_prec)
            for j, (name, start, stop, block) in enumerate(block_spans):
                inv_cov = block_inv_cov[name]
                s2 = max(sigma[j] ** 2, 1e-8)
                if inv_cov is None:
                    prior_precision[start:stop, start:stop] += np.eye(stop - start) / s2
                else:
                    prior_precision[start:stop, start:stop] += inv_cov / s2
            rhs = X.T @ kappa + fixed_prec * fixed_mu
            beta = _chol_solve(precision + prior_precision, rhs, rng)

            # Slice update of each block scale, on the log scale.
            for j, (name, start, stop, block) in enumerate(block_spans):
                b = beta[start:stop]
                inv_cov = block_inv_cov[name]
                quad = float(b @ (inv_cov @ b)) if inv_cov is not None else float(b @ b)
                k_dim = stop - start
                prior_scale = block_prior_scales[j]

                def log_post(s: float, _k=k_dim, _q=quad, _p=prior_scale) -> float:
                    if s <= 0:
                        return -np.inf
                    return (
                        -_k * math.log(s)
                        - 0.5 * _q / (s ** 2)
                        + _half_normal_logpdf(s, _p)
                    )

                sigma[j] = _slice_sample_log(rng, log_post, sigma[j])
                proposals[j] += 1
                accepted[j] += 1

            # Interweaving (ASIS) step: refresh each scale again in the
            # non-centred parameterization, holding the standardized effects
            # z = beta_block / sigma fixed. The centred update above mixes well
            # when a block is strongly identified; this one mixes well when it
            # is weakly identified, which is exactly the funnel that stalls the
            # host-by-family interaction. Interweaving the two is uniformly at
            # least as good as either alone.
            eta_total = X @ beta
            for j, (name, start, stop, block) in enumerate(block_spans):
                s_now = float(sigma[j])
                if s_now <= 1e-10:
                    continue
                z = beta[start:stop] / s_now
                unit = z[block.index]                   # per-observation effect at sigma = 1
                eta_rest = eta_total - s_now * unit
                prior_scale = block_prior_scales[j]

                def log_post_nc(s: float, _rest=eta_rest, _u=unit, _p=prior_scale) -> float:
                    if s <= 0:
                        return -np.inf
                    eta_try = _rest + s * _u
                    ll = float(np.sum(-np.logaddexp(0.0, -np.where(y > 0.5, eta_try, -eta_try))))
                    return ll + _half_normal_logpdf(s, _p)

                s_new = _slice_sample_log(rng, log_post_nc, s_now)
                sigma[j] = s_new
                beta[start:stop] = s_new * z
                eta_total = eta_rest + s_new * unit

            if it >= tune and ((it - tune) % thin == 0) and store < kept:
                beta_draws[store] = beta
                if n_blocks:
                    sigma_draws[store] = sigma
                if store_log_lik:
                    eta_s = X @ beta
                    loglik_draws[store] = -np.logaddexp(0.0, -np.where(y > 0.5, eta_s, -eta_s))
                store += 1

            if progress_every and it % progress_every == 0:
                print(f"[evonexus] chain {chain} iteration {it}/{total}", flush=True)

        chain_beta.append(beta_draws[:store])
        chain_sigma.append(sigma_draws[:store] if n_blocks else np.zeros((store, 0)))
        if store_log_lik:
            chain_loglik.append(loglik_draws[:store])

    posterior = {
        "beta": np.stack(chain_beta),                       # (chain, draw, p)
        "sigma": np.stack(chain_sigma) if n_blocks else None,
        "log_likelihood": np.stack(chain_loglik) if store_log_lik else None,
        "coefficient_names": names,
        "spans": spans,
        "sigma_names": [b[0] for b in block_spans],
        "engine": ENGINE,
        "n_chains": chains,
        "n_draws": chain_beta[0].shape[0],
        "seed": seed,
        "tune": tune,
        "sampler_stats": {
            "divergences": 0,
            "divergence_note": (
                "The Gibbs sampler has no divergence concept. Convergence is "
                "assessed by R-hat and ESS only; treat the absence of "
                "divergences as uninformative, not as evidence of good geometry."
            ),
            "scale_update": "slice sampling on log(sigma); no acceptance rate applies",
            "polya_gamma_truncation": _PG_TRUNCATION,
        },
    }
    return posterior
