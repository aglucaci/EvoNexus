"""Engine selection, fitting and persistence of posterior draws."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from evonexus.bayesian import fallback as fallback_engine
from evonexus.bayesian import model as pymc_engine
from evonexus.bayesian.design import ModelDesign
from evonexus.core.io import read_json, write_json


def available_engines() -> list[str]:
    engines = [fallback_engine.ENGINE]
    if pymc_engine.pymc_available():
        engines.insert(0, pymc_engine.ENGINE)
    return engines


def resolve_engine(requested: str = "auto") -> str:
    """Pick the sampler. ``auto`` prefers PyMC/NUTS when it is installed."""
    if requested in {"auto", None, ""}:
        return available_engines()[0]
    if requested == pymc_engine.ENGINE and not pymc_engine.pymc_available():
        raise RuntimeError(
            "PyMC/ArviZ are not installed but engine='pymc_nuts' was requested. "
            "Install them (see environment.yml) or use engine='auto'."
        )
    if requested not in {pymc_engine.ENGINE, fallback_engine.ENGINE}:
        raise ValueError(f"unknown engine {requested!r}; available: {available_engines()}")
    return requested


def fit_model(
    design: ModelDesign,
    priors: dict[str, Any],
    sampler: dict[str, Any] | None = None,
    engine: str = "auto",
    store_log_lik: bool = True,
) -> dict[str, Any]:
    """Fit one model and return posterior draws plus sampler metadata."""
    sampler = dict(sampler or {})
    chosen = resolve_engine(engine)
    draws = int(sampler.get("draws", 1000))
    tune = int(sampler.get("tune", 1000))
    chains = int(sampler.get("chains", 4))
    seed = int(sampler.get("seed", 20250101))

    if chosen == pymc_engine.ENGINE:
        posterior = pymc_engine.sample_pymc(
            design,
            priors,
            draws=draws,
            tune=tune,
            chains=chains,
            seed=seed,
            target_accept=float(sampler.get("target_accept", 0.9)),
            max_treedepth=int(sampler.get("max_treedepth", 12)),
            cores=sampler.get("cores"),
            store_log_lik=store_log_lik,
        )
    else:
        posterior = fallback_engine.sample(
            design,
            priors,
            draws=draws,
            tune=tune,
            chains=chains,
            seed=seed,
            store_log_lik=store_log_lik,
        )

    posterior["design_summary"] = design.summary()
    posterior["design_metadata"] = design.metadata
    posterior["engine_available"] = available_engines()
    # The same descriptive metadata that save_posterior writes to disk, so an
    # in-memory posterior and a reloaded one behave identically in
    # PosteriorPredictor. Without this, predictions silently drop every fixed
    # effect and treat every group level as unseen.
    posterior.update(_predictor_metadata(design))
    return posterior


def _predictor_metadata(design: ModelDesign) -> dict[str, Any]:
    return {
        "fixed_names": list(design.fixed_names),
        "has_intercept": bool(design.has_intercept),
        "blocks": {
            b.name: {"levels": list(b.levels), "prior_key": b.prior_key}
            for b in design.blocks
        },
    }


# ---------------------------------------------------------------------------
# persistence
# ---------------------------------------------------------------------------
def save_posterior(posterior: dict[str, Any], design: ModelDesign,
                   out_dir: str | Path, name: str) -> dict[str, str]:
    """Persist draws as .npz plus a JSON sidecar describing them."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    arrays = {"beta": np.asarray(posterior["beta"], dtype=np.float32)}
    if posterior.get("sigma") is not None:
        arrays["sigma"] = np.asarray(posterior["sigma"], dtype=np.float32)
    if posterior.get("log_likelihood") is not None:
        arrays["log_likelihood"] = np.asarray(posterior["log_likelihood"], dtype=np.float32)

    npz_path = out / f"posterior_{name}.npz"
    np.savez_compressed(npz_path, **arrays)

    meta = {
        "name": name,
        "engine": posterior.get("engine"),
        "engine_available": posterior.get("engine_available"),
        "coefficient_names": posterior.get("coefficient_names"),
        "spans": [list(s) for s in posterior.get("spans", [])],
        "sigma_names": posterior.get("sigma_names"),
        "n_chains": posterior.get("n_chains"),
        "n_draws": posterior.get("n_draws"),
        "seed": posterior.get("seed"),
        "tune": posterior.get("tune"),
        "sampler_stats": posterior.get("sampler_stats"),
        "design_summary": posterior.get("design_summary"),
        "design_metadata": posterior.get("design_metadata"),
        "blocks": {b.name: {"levels": b.levels, "prior_key": b.prior_key} for b in design.blocks},
        "fixed_names": design.fixed_names,
        "has_intercept": design.has_intercept,
        "arrays": {k: list(v.shape) for k, v in arrays.items()},
    }
    meta_path = write_json(meta, out / f"posterior_{name}.json")

    idata = posterior.get("_idata")
    netcdf_path = None
    if idata is not None:
        try:
            netcdf_path = out / f"idata_{name}.nc"
            idata.to_netcdf(netcdf_path)
        except Exception:
            netcdf_path = None

    result = {"draws": str(npz_path), "metadata": str(meta_path)}
    if netcdf_path:
        result["inference_data"] = str(netcdf_path)
    return result


def load_posterior(out_dir: str | Path, name: str) -> dict[str, Any]:
    out = Path(out_dir)
    meta = read_json(out / f"posterior_{name}.json")
    with np.load(out / f"posterior_{name}.npz") as handle:
        arrays = {k: handle[k].astype(float) for k in handle.files}
    posterior = dict(meta)
    posterior.update(arrays)
    posterior["spans"] = [tuple(s) for s in meta.get("spans", [])]
    return posterior
