# evonexus-bayes

Hierarchical Bayesian modelling of library-level viral RNA detection: the
statistical core extracted from [EvoNexus](../../README.md) so that it can be
installed, versioned and cited on its own.

The package knows nothing about EvoNexus. It takes a `pandas.DataFrame` and an
already-resolved priors mapping, and returns posterior draws, diagnostics,
standardized predictions and contrasts. It reads no configuration files and
owns no output-directory layout — those stay on the pipeline side, in
`scripts/evonexus/bayesian/run.py`.

## Install

```bash
pip install ./packages/evonexus-bayes          # from an EvoNexus checkout
pip install "evonexus-bayes[pymc]"             # with the NUTS engine
```

Without the `pymc` extra the built-in Pólya-Gamma Gibbs sampler is used. Both
engines fit the same `ModelDesign` and are subject to the same R-hat and ESS
reporting; the engine is recorded in every posterior artifact so a fallback fit
can never be mistaken for a NUTS fit.

## Use

```python
from evonexus_bayes import build_design, fit_model, habitat_contrast

design = build_design(frame, terms=["habitat", "host", "tissue"], model_name="M3")
posterior = fit_model(design, priors, sampler=sampler, engine="auto")
contrast = habitat_contrast(design, posterior)
```

## Contents

| Module | Role |
| --- | --- |
| `design` | Builds `ModelDesign`: indices, effect blocks, standardized design matrix |
| `model` | PyMC/NUTS implementation of the detection model |
| `fallback` | Self-contained Pólya-Gamma Gibbs sampler |
| `fit` | Engine selection, fitting, posterior save/load |
| `predict` | Engine-agnostic prediction, integrating over unseen levels |
| `contrast` | Posterior standardization and the habitat estimand |
| `diagnostics` | R-hat, ESS, HDI, PSIS-LOO, model comparison |
| `checks` | Prior-predictive simulation and parameter recovery |

## License

MIT.
