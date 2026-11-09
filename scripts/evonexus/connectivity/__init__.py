"""Posterior cross-host connectivity and its sampling-preserving nulls."""

from evonexus.connectivity.compute import (
    connectivity_matrix,
    host_lineage_probabilities,
    marginal_null,
    permutation_null,
    run_connectivity,
)

__all__ = [
    "host_lineage_probabilities",
    "connectivity_matrix",
    "marginal_null",
    "permutation_null",
    "run_connectivity",
]
