"""Posterior prediction from saved draws, engine-agnostic.

Everything downstream — habitat standardization, connectivity, calibration,
held-out scoring — goes through :class:`PosteriorPredictor`. It reads the flat
coefficient array plus the block/level metadata written by
``evonexus.bayesian.fit.save_posterior``, so a result can be re-derived from
files alone without re-fitting and without knowing which sampler produced them.

Unseen levels (a new sampling site, a new sequencing batch) are *integrated
over* their estimated population distribution rather than estimated from the
data being predicted. That is what makes held-out-site prediction honest.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd


def expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500, 500)))


class PosteriorPredictor:
    """Compute linear predictors for arbitrary covariate profiles."""

    def __init__(self, posterior: Mapping[str, Any]):
        beta = np.asarray(posterior["beta"], dtype=float)
        if beta.ndim == 3:
            self.n_chains, self.n_draws = beta.shape[0], beta.shape[1]
            beta = beta.reshape(-1, beta.shape[-1])
        else:
            self.n_chains, self.n_draws = 1, beta.shape[0]
        self.beta = beta                                  # (S, p)
        self.n_samples = beta.shape[0]

        sigma = posterior.get("sigma")
        self.sigma = (
            np.asarray(sigma, dtype=float).reshape(-1, np.asarray(sigma).shape[-1])
            if sigma is not None and np.asarray(sigma).size else None
        )
        self.sigma_names: list[str] = list(posterior.get("sigma_names") or [])

        self.spans = [tuple(s) for s in posterior.get("spans", [])]
        self.fixed_names: list[str] = list(posterior.get("fixed_names") or [])
        self.has_intercept = bool(posterior.get("has_intercept", True))
        blocks = posterior.get("blocks") or {}
        self.block_levels: dict[str, list[str]] = {
            name: list(spec["levels"]) for name, spec in blocks.items()
        }
        self.block_span: dict[str, tuple[int, int]] = {}
        self.fixed_span: dict[str, tuple[int, int]] = {}
        self.intercept_span: tuple[int, int] | None = None
        for label, start, stop in self.spans:
            if label == "intercept":
                self.intercept_span = (start, stop)
            elif label.startswith("fixed:"):
                self.fixed_span[label.split(":", 1)[1]] = (start, stop)
            elif label.startswith("block:"):
                self.block_span[label.split(":", 1)[1]] = (start, stop)
        self.engine = posterior.get("engine")

    # -- components ------------------------------------------------------
    def intercept(self) -> np.ndarray:
        if not self.has_intercept or self.intercept_span is None:
            return np.zeros(self.n_samples)
        start, _stop = self.intercept_span
        return self.beta[:, start]

    def fixed(self, name: str) -> np.ndarray:
        span = self.fixed_span.get(name)
        if span is None:
            return np.zeros(self.n_samples)
        return self.beta[:, span[0]]

    def block_effects(self, name: str) -> np.ndarray | None:
        span = self.block_span.get(name)
        if span is None:
            return None
        return self.beta[:, span[0]: span[1]]            # (S, n_levels)

    def block_scale(self, name: str) -> np.ndarray:
        """Posterior draws of a block's scale, for integrating over new levels."""
        if self.sigma is not None and name in self.sigma_names:
            return self.sigma[:, self.sigma_names.index(name)]
        effects = self.block_effects(name)
        if effects is None or effects.shape[1] < 2:
            return np.zeros(self.n_samples)
        return effects.std(axis=1, ddof=1)

    def level_index(self, block: str, labels: Sequence[str]) -> np.ndarray:
        levels = self.block_levels.get(block, [])
        lookup = {lvl: i for i, lvl in enumerate(levels)}
        return np.array([lookup.get(str(lbl), -1) for lbl in labels], dtype=int)

    def block_contribution(
        self,
        block: str,
        labels: Sequence[str],
        integrate_unseen: bool = True,
        rng: np.random.Generator | None = None,
    ) -> np.ndarray:
        """Per-draw contribution of a block for a list of level labels.

        Levels absent from the fit are drawn from ``N(0, sigma_block)`` once per
        posterior draw — the correct treatment for an unseen site or batch.
        """
        effects = self.block_effects(block)
        n = len(labels)
        if effects is None:
            return np.zeros((self.n_samples, n))
        idx = self.level_index(block, labels)
        out = np.zeros((self.n_samples, n))
        seen = idx >= 0
        if seen.any():
            out[:, seen] = effects[:, idx[seen]]
        unseen = ~seen
        if unseen.any():
            if not integrate_unseen:
                out[:, unseen] = 0.0
            else:
                rng = rng or np.random.default_rng(0)
                scale = self.block_scale(block)[:, None]
                out[:, unseen] = rng.standard_normal((self.n_samples, int(unseen.sum()))) * scale
        return out

    # -- full linear predictor -------------------------------------------
    def eta(
        self,
        profiles: pd.DataFrame,
        integrate: Iterable[str] = ("site", "batch"),
        rng: np.random.Generator | None = None,
        habitat_override: float | None = None,
    ) -> np.ndarray:
        """Linear predictor of shape (n_samples, n_profiles).

        ``profiles`` carries one column per block (level labels) and per fixed
        effect (numeric). Blocks named in ``integrate`` are always drawn from
        their population distribution, even if the label exists in the fit:
        that is how predictions are standardized to a *new* site.
        """
        rng = rng or np.random.default_rng(20250101)
        n = len(profiles)
        total = np.repeat(self.intercept()[:, None], n, axis=1)
        integrate = set(integrate)

        for name in self.fixed_names:
            if name == "habitat" and habitat_override is not None:
                values = np.full(n, float(habitat_override))
            elif name.startswith("ecology:"):
                column = name.split(":", 1)[1]
                values = pd.to_numeric(profiles.get(column, 0.0), errors="coerce")
                values = np.asarray(pd.Series(values).fillna(0.0), dtype=float) if np.ndim(values) else np.zeros(n)
            else:
                column = {"habitat": "U", "depth": "log_depth_z",
                          "pool_size": "log_pool_size_z", "time": "time_z"}.get(name, name)
                if column not in profiles.columns:
                    values = np.zeros(n)
                else:
                    values = np.asarray(pd.to_numeric(profiles[column], errors="coerce").fillna(0.0), dtype=float)
            total = total + self.fixed(name)[:, None] * values[None, :]

        column_for_block = {
            "virus_family": "virus_family",
            "virus_lineage": "virus_unit",
            "host": "host_species",
            "host_phylo": "host_species",
            "tissue": "tissue_group",
            "site": "site_id",
            "batch": "sequencing_batch",
        }
        for block in self.block_span:
            if block == "host_family_ix":
                labels = (
                    profiles["host_species"].astype(str) + " x " + profiles["virus_family"].astype(str)
                ).tolist()
            else:
                column = column_for_block.get(block, block)
                labels = profiles[column].astype(str).tolist() if column in profiles.columns else ["__absent__"] * n
            force = block in integrate
            if force:
                scale = self.block_scale(block)[:, None]
                contribution = rng.standard_normal((self.n_samples, n)) * scale
            else:
                contribution = self.block_contribution(block, labels, rng=rng)
            total = total + contribution

        return total

    def probability(self, profiles: pd.DataFrame, **kwargs: Any) -> np.ndarray:
        return expit(self.eta(profiles, **kwargs))
