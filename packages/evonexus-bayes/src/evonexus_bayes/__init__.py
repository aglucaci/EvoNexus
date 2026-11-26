"""Bayesian hierarchical modelling of library-level viral RNA detection."""

from importlib.metadata import PackageNotFoundError, version

from evonexus_bayes.contrast import habitat_contrast, standardized_predictions
from evonexus_bayes.design import ModelDesign, build_design
from evonexus_bayes.diagnostics import convergence_report, summarize_posterior
from evonexus_bayes.fit import available_engines, fit_model

try:
    __version__ = version("evonexus-bayes")
except PackageNotFoundError:  # running from the source tree, not installed
    __version__ = "0.1.0"

__all__ = [
    "ModelDesign",
    "build_design",
    "fit_model",
    "available_engines",
    "summarize_posterior",
    "convergence_report",
    "habitat_contrast",
    "standardized_predictions",
    "__version__",
]
