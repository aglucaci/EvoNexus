"""Resolve analysis artifacts and report honestly on what is missing.

Nothing in the reporting layer invents a value. If an artifact is absent, the
figure is emitted as a labelled placeholder, the table is omitted with a
reason, and the legend says "not generated" instead of a sample size.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from evonexus.core.io import read_json, read_table, resolve_table


@dataclass
class ArtifactStore:
    """Lazy, forgiving access to the results tree."""

    root: Path
    missing: List[str] = field(default_factory=list)
    _cache: Dict[str, Any] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.root = Path(self.root)

    # -- low-level -------------------------------------------------------
    def path(self, *parts: str) -> Path:
        return self.root.joinpath(*parts)

    def table(self, *parts: str) -> Optional[pd.DataFrame]:
        key = "table:" + "/".join(parts)
        if key in self._cache:
            return self._cache[key]
        try:
            df = read_table(self.path(*parts))
        except (FileNotFoundError, OSError):
            self.missing.append("/".join(parts))
            df = None
        self._cache[key] = df
        return df

    def json(self, *parts: str) -> Optional[Dict[str, Any]]:
        key = "json:" + "/".join(parts)
        if key in self._cache:
            return self._cache[key]
        p = self.path(*parts)
        if not p.exists():
            self.missing.append("/".join(parts))
            doc = None
        else:
            try:
                doc = read_json(p)
            except Exception:
                self.missing.append("/".join(parts))
                doc = None
        self._cache[key] = doc
        return doc

    def exists(self, *parts: str) -> bool:
        p = self.path(*parts)
        if p.exists():
            return True
        try:
            resolve_table(p)
            return True
        except FileNotFoundError:
            return False

    # -- named artifacts --------------------------------------------------
    def canonical(self, name: str) -> Optional[pd.DataFrame]:
        return self.table("processed", name)

    def qc(self) -> Optional[Dict[str, Any]]:
        return self.json("qc", "qc_gates.json")

    def waterfall(self) -> Optional[pd.DataFrame]:
        return self.table("qc", "exclusion_waterfall")

    def analysis_frame(self, resolution: str = "lineage") -> Optional[pd.DataFrame]:
        return self.table("matrix", f"analysis_frame_{resolution}")

    def matrix_metadata(self, resolution: str = "lineage") -> Optional[Dict[str, Any]]:
        return self.json("matrix", f"matrix_metadata_{resolution}.json")

    def contrast(self) -> Optional[pd.DataFrame]:
        return self.table("bayesian", "habitat_contrast")

    def contrast_by_family(self) -> Optional[pd.DataFrame]:
        return self.table("bayesian", "habitat_contrast_by_family")

    def contrast_draws(self) -> Optional[pd.DataFrame]:
        return self.table("bayesian", "habitat_contrast_draws")

    def contrast_metadata(self) -> Optional[Dict[str, Any]]:
        return self.json("bayesian", "habitat_contrast_metadata.json")

    def raw_rates(self) -> Optional[pd.DataFrame]:
        return self.table("bayesian", "raw_detection_rates")

    def sensitivity(self) -> Optional[pd.DataFrame]:
        return self.table("bayesian", "sensitivity")

    def model_comparison(self) -> Optional[pd.DataFrame]:
        return self.table("bayesian", "model_comparison")

    def bayesian_status(self) -> Optional[Dict[str, Any]]:
        return self.json("bayesian", "bayesian_status.json")

    def prior_check(self) -> Optional[Dict[str, Any]]:
        return self.json("bayesian", "prior_predictive_check.json")

    def posterior_predictive(self) -> Optional[Dict[str, Any]]:
        return self.json("bayesian", "posterior_predictive_check.json")

    def parameter_recovery(self) -> Optional[pd.DataFrame]:
        return self.table("bayesian", "parameter_recovery")

    def convergence(self, model: str) -> Optional[Dict[str, Any]]:
        return self.json("bayesian", f"convergence_{model}.json")

    def posterior_summary(self, model: str) -> Optional[pd.DataFrame]:
        return self.table("bayesian", f"posterior_summary_{model}")

    def connectivity(self) -> Optional[pd.DataFrame]:
        return self.table("connectivity", "connectivity")

    def connectivity_metadata(self) -> Optional[Dict[str, Any]]:
        return self.json("connectivity", "connectivity_metadata.json")

    def connectivity_permutation(self) -> Optional[pd.DataFrame]:
        return self.table("connectivity", "connectivity_permutation_null")

    def selection_results(self) -> Optional[pd.DataFrame]:
        return self.table("evolution", "selection_results")

    def evolution_status(self) -> Optional[Dict[str, Any]]:
        return self.json("evolution", "evolution_status.json")

    def lineage_results(self) -> Optional[pd.DataFrame]:
        return self.table("evolution", "evolution_lineage_results")

    def cv_folds(self) -> Optional[pd.DataFrame]:
        return self.table("validation", "cv_folds")

    def calibration(self, group: str = "site_id") -> Optional[pd.DataFrame]:
        return self.table("validation", f"calibration_{group}")

    def validation_metadata(self) -> Optional[Dict[str, Any]]:
        return self.json("validation", "validation_metadata.json")

    def manifest(self) -> Optional[Dict[str, Any]]:
        return self.json("acquisition", "manifest.json")

    # -- provenance -------------------------------------------------------
    def is_synthetic(self) -> bool:
        for doc in (self.bayesian_status(), self.qc(), self.matrix_metadata(),
                    self.manifest()):
            if isinstance(doc, dict) and doc.get("synthetic"):
                return True
        return False

    def availability(self, requirements: Dict[str, bool]) -> Dict[str, Any]:
        present = {k: v for k, v in requirements.items() if v}
        absent = sorted(k for k, v in requirements.items() if not v)
        return {"available": sorted(present), "missing": absent, "complete": not absent}
