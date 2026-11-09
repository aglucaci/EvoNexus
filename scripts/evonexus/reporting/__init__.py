"""Manuscript outputs: figures, source tables and generated legends."""

from evonexus.reporting.artifacts import ArtifactStore
from evonexus.reporting.figures import build_figures
from evonexus.reporting.legends import build_legends
from evonexus.reporting.run import run_reporting
from evonexus.reporting.tables import build_tables

__all__ = [
    "ArtifactStore",
    "build_figures",
    "build_tables",
    "build_legends",
    "run_reporting",
]
