"""Library behaviour, checked without the pipeline.

Everything here constructs its own data and passes priors as literals. Tests
that assert the *repository's* shipped `config/priors.yaml` drives the model,
or that depend on the harmonized analysis frame, are integration tests and
stay in the EvoNexus suite.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from evonexus_bayes.design import build_design
from evonexus_bayes.diagnostics import ess, hdi, rhat
from evonexus_bayes.fit import fit_model
from evonexus_bayes.predict import PosteriorPredictor, expit


def test_rhat_and_ess_behave_on_known_input():
    rng = np.random.default_rng(0)
    independent = rng.standard_normal((4, 500))
    assert rhat(independent) < 1.05
    assert ess(independent) > 400
    # Chains stuck at different values must be flagged.
    stuck = np.stack([np.full(500, -3.0), np.full(500, 3.0)])
    assert rhat(stuck) > 1.5 or np.isnan(rhat(stuck))


def test_hdi_covers_the_bulk():
    rng = np.random.default_rng(1)
    draws = rng.normal(2.0, 1.0, 20000)
    low, high = hdi(draws, 0.95)
    assert low < 2.0 < high
    assert 3.5 < (high - low) < 4.5


def test_recovers_a_known_logistic_coefficient():
    """End-to-end check that the sampler is not biased on a simple case."""
    rng = np.random.default_rng(5)
    n = 4000
    x = rng.standard_normal(n)
    group = rng.integers(0, 6, n)
    offsets = rng.normal(0, 0.5, 6)
    eta = -1.0 + 1.25 * x + offsets[group]
    y = (rng.random(n) < expit(eta)).astype(float)

    frame = pd.DataFrame(
        {
            "y": y,
            "U": x,                      # reuse the habitat slot as the known covariate
            "virus_unit": "V",
            "virus_family": "F",
            "host_species": "H",
            "tissue_group": "T",
            "site_id": [f"S{g}" for g in group],
            "sequencing_batch": "B",
            "log_depth_z": 0.0,
            "log_pool_size_z": 0.0,
            "time_z": 0.0,
        }
    )
    design = build_design(frame, ["intercept", "habitat", "site"], model_name="known")
    posterior = fit_model(design, {"intercept": {"dist": "normal", "mu": 0, "sigma": 2},
                                   "habitat_beta": {"dist": "normal", "mu": 0, "sigma": 2},
                                   "site_sd": {"dist": "halfnormal", "sigma": 1}},
                          sampler={"draws": 400, "tune": 400, "chains": 2, "seed": 9},
                          engine="polya_gamma_gibbs")
    predictor = PosteriorPredictor(posterior)
    beta = predictor.fixed("habitat")
    low, high = hdi(beta, 0.95)
    assert low <= 1.25 <= high, f"true 1.25 outside HDI [{low:.3f}, {high:.3f}]"
    assert abs(float(np.mean(beta)) - 1.25) < 0.25
