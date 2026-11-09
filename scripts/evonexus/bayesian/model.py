"""PyMC implementation of the hierarchical detection model.

Likelihood
----------
For library :math:`l` and viral lineage unit :math:`v`,

.. math::

    Y_{lv} \\sim \\mathrm{Bernoulli}(p_{lv})

.. math::

    \\mathrm{logit}(p_{lv}) = \\alpha + a_{f(v)} + a_v + \\beta_U U_l
        + b_{h(l)} + \\kappa_{h(l),f(v)} + c_{t(l)} + g_{r(l)}
        + \\gamma_D \\log D_l + \\gamma_N \\log N_l + \\gamma_T T_l
        + \\theta^\\top E_l

Only cells that were actually assessed enter the likelihood. Varying effects
use a non-centred parameterization, which is what makes the host-by-family
interaction tractable when most host x family combinations are sparse.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from evonexus.bayesian.design import ModelDesign, design_matrix

ENGINE = "pymc_nuts"


def pymc_available() -> bool:
    try:
        import arviz  # noqa: F401
        import pymc  # noqa: F401
    except Exception:
        return False
    return True


def _dist(pm, name: str, spec: dict[str, Any], shape=None):
    kind = str(spec.get("dist", "normal")).lower()
    if kind == "normal":
        return pm.Normal(name, mu=float(spec.get("mu", 0.0)),
                         sigma=float(spec.get("sigma", 1.0)), shape=shape)
    if kind in {"halfnormal", "half_normal"}:
        return pm.HalfNormal(name, sigma=float(spec.get("sigma", 1.0)), shape=shape)
    if kind == "beta":
        return pm.Beta(name, alpha=float(spec.get("alpha", 2.0)),
                       beta=float(spec.get("beta", 2.0)), shape=shape)
    if kind in {"exponential", "exp"}:
        return pm.Exponential(name, lam=float(spec.get("lam", 1.0)), shape=shape)
    raise ValueError(f"unsupported prior distribution {kind!r}")


def build_pymc_model(design: ModelDesign, priors: dict[str, Any]):
    """Construct the PyMC model object for a design."""
    import pymc as pm
    import pytensor.tensor as pt

    coords = {name: levels for name, levels in design.coords.items()}
    coords["obs"] = np.arange(design.n_obs)
    if design.fixed_names:
        coords["fixed"] = design.fixed_names

    with pm.Model(coords=coords) as model:
        eta = 0.0

        if design.has_intercept:
            alpha = _dist(pm, "intercept", priors.get("intercept", {}))
            eta = eta + alpha

        if design.n_fixed:
            # Each fixed effect keeps its own prior rather than sharing one.
            parts = []
            for i, name in enumerate(design.fixed_names):
                spec = priors.get(design.fixed_prior_keys[i], {})
                safe = name.replace(":", "_")
                parts.append(_dist(pm, f"b_{safe}", spec))
            theta = pt.stack(parts)
            pm.Deterministic("fixed_effects", theta, dims="fixed")
            eta = eta + pt.dot(pm.math.constant(design.X_fixed), theta)

        for block in design.blocks:
            sd = _dist(pm, f"sigma_{block.name}", priors.get(block.prior_key, {}))
            offset = pm.Normal(f"z_{block.name}", 0.0, 1.0, dims=block.name)
            if block.covariance is not None and not np.allclose(block.covariance, np.eye(block.size)):
                lam_spec = priors.get("host_phylo_lambda", {"dist": "beta", "alpha": 2.0, "beta": 2.0})
                lam = _dist(pm, f"lambda_{block.name}", lam_spec)
                mixed = lam * block.covariance + (1.0 - lam) * np.eye(block.size)
                chol = pt.slinalg.cholesky(mixed + 1e-6 * np.eye(block.size))
                effect = pm.Deterministic(block.name, sd * pt.dot(chol, offset), dims=block.name)
            else:
                effect = pm.Deterministic(block.name, sd * offset, dims=block.name)
            eta = eta + effect[block.index]

        p = pm.Deterministic("p", pm.math.sigmoid(eta), dims="obs")
        pm.Bernoulli("y", p=p, observed=design.y, dims="obs")
    return model


def _flatten_posterior(idata, design: ModelDesign, names: list[str]) -> np.ndarray:
    """Stack PyMC posterior variables into the design-matrix column order."""
    post = idata.posterior
    pieces: list[np.ndarray] = []

    if design.has_intercept:
        pieces.append(np.asarray(post["intercept"])[..., None])

    for i, name in enumerate(design.fixed_names):
        safe = name.replace(":", "_")
        pieces.append(np.asarray(post[f"b_{safe}"])[..., None])

    for block in design.blocks:
        pieces.append(np.asarray(post[block.name]))

    return np.concatenate(pieces, axis=-1)


def sample_pymc(
    design: ModelDesign,
    priors: dict[str, Any],
    draws: int = 1000,
    tune: int = 1000,
    chains: int = 4,
    seed: int = 20250101,
    target_accept: float = 0.9,
    max_treedepth: int = 12,
    cores: int | None = None,
    store_log_lik: bool = True,
) -> dict[str, Any]:
    """Fit with NUTS and return the same dict structure as the fallback."""
    import pymc as pm

    _X, names, spans = design_matrix(design)
    model = build_pymc_model(design, priors)
    with model:
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            random_seed=seed,
            target_accept=target_accept,
            max_treedepth=max_treedepth,
            cores=cores if cores is not None else min(chains, 4),
            idata_kwargs={"log_likelihood": store_log_lik},
            progressbar=False,
        )

    beta = _flatten_posterior(idata, design, names)
    sigma_names = [b.name for b in design.blocks]
    sigma = (
        np.stack([np.asarray(idata.posterior[f"sigma_{n}"]) for n in sigma_names], axis=-1)
        if sigma_names else None
    )
    log_lik = None
    if store_log_lik and hasattr(idata, "log_likelihood"):
        log_lik = np.asarray(idata.log_likelihood["y"])

    stats = getattr(idata, "sample_stats", None)
    divergences = int(np.asarray(stats["diverging"]).sum()) if stats is not None and "diverging" in stats else 0
    treedepth_hits = 0
    if stats is not None and "tree_depth" in stats:
        treedepth_hits = int((np.asarray(stats["tree_depth"]) >= max_treedepth).sum())

    return {
        "beta": beta,
        "sigma": sigma,
        "log_likelihood": log_lik,
        "coefficient_names": names,
        "spans": spans,
        "sigma_names": sigma_names,
        "engine": ENGINE,
        "n_chains": chains,
        "n_draws": beta.shape[1],
        "seed": seed,
        "tune": tune,
        "sampler_stats": {
            "divergences": divergences,
            "divergence_fraction": divergences / max(1, beta.shape[0] * beta.shape[1]),
            "treedepth_hits": treedepth_hits,
            "target_accept": target_accept,
            "max_treedepth": max_treedepth,
        },
        "_idata": idata,
    }


def prior_predictive(design: ModelDesign, priors: dict[str, Any], draws: int = 1000,
                     seed: int = 20250101) -> np.ndarray:
    """Simulate detections under the priors alone (PyMC path)."""
    import pymc as pm

    model = build_pymc_model(design, priors)
    with model:
        pp = pm.sample_prior_predictive(draws=draws, random_seed=seed)
    return np.asarray(pp.prior_predictive["y"])
