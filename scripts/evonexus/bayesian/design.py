"""Model design: indices, blocks and the standardized design matrix.

The design is engine-agnostic. Both the PyMC path and the self-contained
Gibbs fallback consume the same :class:`ModelDesign`, so a model specified
once is fitted identically whichever sampler is available, and the posterior
standardization code does not need to know which one ran.

Effect blocks
-------------
Each varying-effect term becomes a *block* of indicator columns sharing one
scale parameter. Fixed effects are single columns with a prior taken straight
from ``config/priors.yaml``. The block structure is what makes the
host-by-family interaction shrinkable rather than a free parameter per cell.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

# term -> (kind, prior key)
TERM_SPEC: Dict[str, Tuple[str, str]] = {
    "intercept": ("intercept", "intercept"),
    "virus_family": ("block", "family_sd"),
    "virus_lineage": ("block", "lineage_sd"),
    "host": ("block", "host_sd"),
    "tissue": ("block", "tissue_sd"),
    "site": ("block", "site_sd"),
    "batch": ("block", "batch_sd"),
    "host_family_ix": ("block", "interaction_sd"),
    "habitat": ("fixed", "habitat_beta"),
    "depth": ("fixed", "depth_gamma"),
    "pool_size": ("fixed", "poolsize_gamma"),
    "time": ("fixed", "time_gamma"),
    "ecology": ("fixed_multi", "ecology_theta"),
    "host_phylo": ("block_cov", "host_phylo_sd"),
}

TERM_COLUMN: Dict[str, str] = {
    "virus_family": "virus_family",
    "virus_lineage": "virus_unit",
    "host": "host_species",
    "tissue": "tissue_group",
    "site": "site_id",
    "batch": "sequencing_batch",
    "habitat": "U",
    "depth": "log_depth_z",
    "pool_size": "log_pool_size_z",
    "time": "time_z",
}


@dataclass
class Block:
    """A group of indicator columns sharing one hierarchical scale."""

    name: str
    prior_key: str
    levels: List[str]
    index: np.ndarray                 # per-observation level index
    covariance: np.ndarray | None = None   # for phylogenetic blocks

    @property
    def size(self) -> int:
        return len(self.levels)


@dataclass
class ModelDesign:
    model_name: str
    terms: List[str]
    y: np.ndarray
    n_obs: int
    blocks: List[Block] = field(default_factory=list)
    fixed_names: List[str] = field(default_factory=list)
    fixed_prior_keys: List[str] = field(default_factory=list)
    X_fixed: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    has_intercept: bool = True
    coords: Dict[str, List[str]] = field(default_factory=dict)
    frame: pd.DataFrame | None = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_fixed(self) -> int:
        return self.X_fixed.shape[1] if self.X_fixed.size else 0

    def block(self, name: str) -> Block | None:
        for b in self.blocks:
            if b.name == name:
                return b
        return None

    def parameter_count(self) -> int:
        return int(self.has_intercept) + self.n_fixed + sum(b.size for b in self.blocks)

    def summary(self) -> Dict[str, Any]:
        return {
            "model": self.model_name,
            "terms": self.terms,
            "n_obs": self.n_obs,
            "n_parameters": self.parameter_count(),
            "blocks": {b.name: b.size for b in self.blocks},
            "fixed_effects": self.fixed_names,
            "detection_rate": float(np.mean(self.y)) if self.n_obs else None,
        }


def _factorize(values: pd.Series) -> Tuple[np.ndarray, List[str]]:
    cat = pd.Categorical(values.astype(str))
    return np.asarray(cat.codes, dtype=np.int64), [str(c) for c in cat.categories]


def build_design(
    frame: pd.DataFrame,
    terms: Sequence[str],
    model_name: str = "model",
    ecology_columns: Sequence[str] | None = None,
    host_covariance: Dict[str, Any] | None = None,
) -> ModelDesign:
    """Build a :class:`ModelDesign` for one member of the model ladder."""
    if frame.empty:
        raise ValueError("cannot build a design from an empty analysis frame")

    terms = [t for t in terms]
    unknown = [t for t in terms if t not in TERM_SPEC]
    if unknown:
        raise KeyError(f"unknown model terms: {unknown}")

    y = pd.to_numeric(frame["y"], errors="coerce").to_numpy(dtype=float)
    if np.isnan(y).any():
        raise ValueError(
            "the analysis frame contains missing responses; unassessable cells "
            "must be removed before the design is built, not imputed"
        )

    blocks: List[Block] = []
    coords: Dict[str, List[str]] = {}
    fixed_cols: List[np.ndarray] = []
    fixed_names: List[str] = []
    fixed_keys: List[str] = []

    ecology_columns = list(ecology_columns or [c for c in frame.columns if c.endswith("_z")
                                               and c not in {"log_depth_z", "log_pool_size_z", "time_z"}])

    # Factorize the grouping variables that several terms share.
    family_idx, family_levels = _factorize(frame["virus_family"])
    host_idx, host_levels = _factorize(frame["host_species"])

    for term in terms:
        kind, prior_key = TERM_SPEC[term]

        if kind == "intercept":
            continue

        if kind == "block":
            if term == "host_family_ix":
                combo = (
                    frame["host_species"].astype(str) + " x " + frame["virus_family"].astype(str)
                )
                idx, levels = _factorize(combo)
            elif term == "virus_family":
                idx, levels = family_idx, family_levels
            elif term == "host":
                idx, levels = host_idx, host_levels
            else:
                column = TERM_COLUMN[term]
                if column not in frame.columns:
                    raise KeyError(f"term {term!r} needs column {column!r}, which is absent")
                idx, levels = _factorize(frame[column])
            if len(levels) < 2:
                # A single-level grouping is unidentified; drop it and record why.
                continue
            blocks.append(Block(term, prior_key, levels, idx))
            coords[term] = levels
            continue

        if kind == "block_cov":
            cov = None
            levels = host_levels
            idx = host_idx
            if host_covariance:
                order = [str(s) for s in host_covariance.get("species", [])]
                matrix = np.asarray(host_covariance.get("matrix"), dtype=float)
                lookup = {s: i for i, s in enumerate(order)}
                if all(s in lookup for s in levels) and matrix.size:
                    sel = [lookup[s] for s in levels]
                    cov = matrix[np.ix_(sel, sel)]
            if cov is None:
                # No usable phylogeny: fall back to an exchangeable host block
                # and say so, rather than silently inventing a covariance.
                cov = np.eye(len(levels))
            blocks.append(Block(term, prior_key, levels, idx, covariance=cov))
            coords[term] = levels
            continue

        if kind == "fixed":
            column = TERM_COLUMN[term]
            if column not in frame.columns:
                raise KeyError(f"term {term!r} needs column {column!r}, which is absent")
            values = pd.to_numeric(frame[column], errors="coerce").fillna(0.0).to_numpy(float)
            if np.allclose(values, values[0]):
                continue
            fixed_cols.append(values)
            fixed_names.append(term)
            fixed_keys.append(prior_key)
            continue

        if kind == "fixed_multi":
            for col in ecology_columns:
                if col not in frame.columns:
                    continue
                values = pd.to_numeric(frame[col], errors="coerce").fillna(0.0).to_numpy(float)
                if np.allclose(values, values[0]):
                    continue
                fixed_cols.append(values)
                fixed_names.append(f"ecology:{col}")
                fixed_keys.append(prior_key)

    X_fixed = np.column_stack(fixed_cols) if fixed_cols else np.zeros((len(y), 0))

    # Drop a redundant host block when the phylogenetic block replaces it.
    if any(b.name == "host_phylo" for b in blocks):
        blocks = [b for b in blocks if b.name != "host"]
        coords.pop("host", None)

    design = ModelDesign(
        model_name=model_name,
        terms=terms,
        y=y,
        n_obs=len(y),
        blocks=blocks,
        fixed_names=fixed_names,
        fixed_prior_keys=fixed_keys,
        X_fixed=X_fixed,
        has_intercept="intercept" in terms,
        coords=coords,
        frame=frame.reset_index(drop=True),
    )
    design.metadata = {
        "ecology_columns_used": [n.split(":", 1)[1] for n in fixed_names if n.startswith("ecology:")],
        "terms_requested": list(terms),
        "terms_dropped": [
            t for t in terms
            if t not in {"intercept"}
            and TERM_SPEC[t][0] in {"block", "block_cov"}
            and design.block(t) is None
            and not (t == "host" and design.block("host_phylo"))
        ],
        "host_phylogeny_used": bool(
            host_covariance and design.block("host_phylo") is not None
            and not np.allclose(design.block("host_phylo").covariance, np.eye(design.block("host_phylo").size))
        ),
    }
    return design


def design_matrix(design: ModelDesign) -> Tuple[np.ndarray, List[str], List[Tuple[str, int, int]]]:
    """Dense design matrix plus column names and (block, start, stop) spans."""
    parts: List[np.ndarray] = []
    names: List[str] = []
    spans: List[Tuple[str, int, int]] = []
    cursor = 0

    if design.has_intercept:
        parts.append(np.ones((design.n_obs, 1)))
        names.append("intercept")
        spans.append(("intercept", cursor, cursor + 1))
        cursor += 1

    if design.n_fixed:
        parts.append(design.X_fixed)
        names.extend(design.fixed_names)
        for i, name in enumerate(design.fixed_names):
            spans.append((f"fixed:{name}", cursor + i, cursor + i + 1))
        cursor += design.n_fixed

    for block in design.blocks:
        indicators = np.zeros((design.n_obs, block.size))
        indicators[np.arange(design.n_obs), block.index] = 1.0
        parts.append(indicators)
        names.extend([f"{block.name}[{lvl}]" for lvl in block.levels])
        spans.append((f"block:{block.name}", cursor, cursor + block.size))
        cursor += block.size

    X = np.hstack(parts) if parts else np.zeros((design.n_obs, 0))
    return X, names, spans
