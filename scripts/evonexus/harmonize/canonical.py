"""Construct the canonical tables from mapped source data.

Design rules enforced here:

* A value that the source does not state becomes an explicit ``unknown``
  category, never a modelled level and never a zero.
* A library x lineage cell is a *negative* only when the library was screened
  against the common reference panel at adequate depth. Otherwise it is
  ``unassessable`` and carries no numeric response.
* A mixed-species pool keeps its pooled identity. It is flagged, not split
  into one animal per constituent species.
* Every transformation writes a crosswalk row, so a canonical identifier can
  be traced back to the source table and column it came from.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from evonexus.core.io import read_any_table
from evonexus.harmonize.mapping import ResolvedMapping, _normalize

UNKNOWN = "unknown"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _to_missing(series: pd.Series, tokens: List[str]) -> pd.Series:
    token_set = {str(t).strip().lower() for t in tokens}
    out = series.copy()
    if out.dtype == object or pd.api.types.is_string_dtype(out):
        mask = out.astype(str).str.strip().str.lower().isin(token_set)
        out = out.mask(mask)
    return out


def _load_entity(mapping: ResolvedMapping, entity: str,
                 raw_dir: Path) -> Tuple[pd.DataFrame | None, pd.DataFrame]:
    """Return (mapped dataframe, crosswalk rows) for one canonical entity."""
    res = mapping.entities.get(entity)
    crosswalk = pd.DataFrame(
        columns=["entity", "source_id", "canonical_id", "source_table",
                 "source_column", "evidence", "provenance"]
    )
    if res is None or res.source_path is None:
        return None, crosswalk

    # The inventory may record an absolute path, a path relative to the
    # current directory, or a path relative to the raw-data root. Try them in
    # that order rather than assuming one of them.
    candidates = [Path(res.source_path)]
    if res.source_table:
        candidates.append(Path(raw_dir) / res.source_table)
    path = next((c for c in candidates if c.exists()), None)
    if path is None:
        raise FileNotFoundError(
            f"source table for {entity!r} not found; tried "
            + ", ".join(str(c) for c in candidates)
        )
    raw = read_any_table(path)
    colmap = res.column_map()
    kept = raw[[c for c in colmap if c in raw.columns]].rename(columns=colmap)
    for col in kept.columns:
        kept[col] = _to_missing(kept[col], mapping.missing_tokens)
    kept.attrs["source_table"] = res.source_table
    kept.attrs["source_path"] = str(path)
    kept.attrs["raw"] = raw
    kept.attrs["colmap"] = colmap
    return kept, crosswalk


def _vocab_map(value: Any, vocab: Dict[str, List[str]], default: str = UNKNOWN) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    text = str(value).strip().lower()
    if not text:
        return default
    for level, synonyms in vocab.items():
        for syn in synonyms:
            if text == str(syn).strip().lower():
                return level
    for level, synonyms in vocab.items():
        for syn in synonyms:
            s = str(syn).strip().lower()
            if s and s in text:
                return level
    return default


def _as_bool(value: Any, vocabularies: Dict[str, Any]) -> bool | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip().lower()
    truthy = {str(v).strip().lower() for v in vocabularies.get("boolean_true", [])}
    falsy = {str(v).strip().lower() for v in vocabularies.get("boolean_false", [])}
    if text in truthy:
        return True
    if text in falsy:
        return False
    return None


def _crosswalk_rows(entity: str, ids: pd.Series, source_table: str | None,
                    source_column: str | None, evidence: str) -> pd.DataFrame:
    values = ids.dropna().astype(str).unique()
    return pd.DataFrame(
        {
            "entity": entity,
            "source_id": values,
            "canonical_id": values,
            "source_table": source_table,
            "source_column": source_column,
            "evidence": evidence,
            "provenance": "evonexus.harmonize.canonical",
        }
    )


# ---------------------------------------------------------------------------
# entity builders
# ---------------------------------------------------------------------------
def _build_sites(df: pd.DataFrame | None, mapping: ResolvedMapping) -> pd.DataFrame:
    cols = ["site_id", "region", "country", "habitat", "habitat_raw",
            "latitude", "longitude", "geo_resolution"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)
    out = pd.DataFrame(index=df.index)
    out["site_id"] = df.get("site_id").astype(str)
    out["region"] = df.get("region")
    out["country"] = df.get("country")
    habitat_vocab = (mapping.vocabularies.get("habitat") or {})
    raw_habitat = df.get("habitat_raw")
    out["habitat_raw"] = raw_habitat
    out["habitat"] = [
        _vocab_map(v, habitat_vocab) for v in (raw_habitat if raw_habitat is not None else [])
    ] if raw_habitat is not None else UNKNOWN
    for c in ("latitude", "longitude"):
        out[c] = pd.to_numeric(df.get(c), errors="coerce") if c in df else np.nan
    out["geo_resolution"] = df.get("geo_resolution")
    out = out.drop_duplicates(subset=["site_id"]).reset_index(drop=True)
    return out[cols]


def _build_hosts(df: pd.DataFrame | None, mapping: ResolvedMapping) -> pd.DataFrame:
    cols = ["host_id", "site_id", "host_species", "host_species_resolved",
            "host_genus", "host_family", "pool_status", "pool_size",
            "mixed_species_pool", "collection_date"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)
    out = pd.DataFrame(index=df.index)
    out["host_id"] = df.get("host_id").astype(str)
    out["site_id"] = df.get("site_id").astype(str) if "site_id" in df else None
    species = df.get("host_species")
    out["host_species"] = species
    resolved = []
    mixed_flags = []
    for v in (species if species is not None else []):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            resolved.append(UNKNOWN)
            mixed_flags.append(False)
            continue
        text = str(v).strip()
        is_mixed = bool(re.search(r"(?i)\b(mixed|multiple|several|spp\.?$|/)\b", text)) or "," in text
        mixed_flags.append(is_mixed)
        if is_mixed or not text:
            resolved.append("mixed" if is_mixed else UNKNOWN)
        else:
            resolved.append(re.sub(r"\s+", " ", text))
    out["host_species_resolved"] = resolved if species is not None else UNKNOWN
    out["host_genus"] = df.get("host_genus")
    if "host_genus" not in df and species is not None:
        out["host_genus"] = [
            s.split(" ")[0] if s not in {UNKNOWN, "mixed"} else UNKNOWN for s in out["host_species_resolved"]
        ]
    out["host_family"] = df.get("host_family")

    pool_size = pd.to_numeric(df.get("pool_size"), errors="coerce") if "pool_size" in df else pd.Series(np.nan, index=df.index)
    explicit_pool = [_as_bool(v, mapping.vocabularies) for v in df.get("is_pool", pd.Series([None] * len(df), index=df.index))]
    status = []
    for flag, n in zip(explicit_pool, pool_size):
        if flag is True or (pd.notna(n) and n > 1):
            status.append("pool")
        elif flag is False or (pd.notna(n) and n == 1):
            status.append("individual")
        else:
            status.append(UNKNOWN)
    out["pool_status"] = status
    out["pool_size"] = pool_size
    out["mixed_species_pool"] = mixed_flags if species is not None else False
    out["collection_date"] = df.get("collection_date")
    out = out.drop_duplicates(subset=["host_id"]).reset_index(drop=True)
    return out[cols]


def _effective_depth(raw: pd.DataFrame, read_depth: pd.Series,
                     filters: Dict[str, Any]) -> Tuple[pd.Series, str]:
    """Pick one depth denominator and use it for every library."""
    spec = filters.get("effective_depth", {}) or {}
    primary = str(spec.get("definition", "non_host_non_rrna_reads"))
    fallback = str(spec.get("fallback_definition", "total_clean_reads"))
    norm_cols = {_normalize(c): c for c in raw.columns}
    for name, label in ((primary, primary), (fallback, fallback)):
        col = norm_cols.get(_normalize(name))
        if col is not None:
            values = pd.to_numeric(raw[col], errors="coerce")
            if values.notna().any():
                return values, label
    return read_depth.astype(float), "read_depth_fallback"


def _build_libraries(df: pd.DataFrame | None, mapping: ResolvedMapping,
                     filters: Dict[str, Any]) -> pd.DataFrame:
    cols = ["library_id", "host_id", "site_id", "tissue", "tissue_group",
            "read_depth", "effective_depth", "effective_depth_definition",
            "sequencing_batch", "sra_run", "bioproject", "collection_date",
            "collection_time_numeric", "screened_panel", "is_negative_control"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)
    raw = df.attrs.get("raw", df)
    out = pd.DataFrame(index=df.index)
    out["library_id"] = df.get("library_id").astype(str)
    out["host_id"] = df.get("host_id").astype(str) if "host_id" in df else None
    out["site_id"] = df.get("site_id").astype(str) if "site_id" in df else None
    tissue = df.get("tissue")
    out["tissue"] = tissue
    tissue_vocab = mapping.vocabularies.get("tissue") or {}
    out["tissue_group"] = [_vocab_map(v, tissue_vocab) for v in (tissue if tissue is not None else [])] if tissue is not None else UNKNOWN

    read_depth = pd.to_numeric(df.get("read_depth"), errors="coerce") if "read_depth" in df else pd.Series(np.nan, index=df.index)
    out["read_depth"] = read_depth
    eff, definition = _effective_depth(raw, read_depth, filters)
    out["effective_depth"] = eff.reindex(df.index) if len(eff) == len(df) else eff
    out["effective_depth_definition"] = definition

    out["sequencing_batch"] = df.get("sequencing_batch")
    out["sra_run"] = df.get("sra_run")
    out["bioproject"] = df.get("bioproject")
    date = df.get("collection_date")
    out["collection_date"] = date
    if date is not None:
        parsed = pd.to_datetime(date, errors="coerce")
        if parsed.notna().any():
            origin = parsed.min()
            out["collection_time_numeric"] = (parsed - origin).dt.days / 365.25
        else:
            out["collection_time_numeric"] = np.nan
    else:
        out["collection_time_numeric"] = np.nan

    panel_col = {_normalize(c): c for c in raw.columns}.get(_normalize("reference_panel"))
    out["screened_panel"] = raw[panel_col].astype(str) if panel_col else UNKNOWN

    pattern = (filters.get("contamination", {}) or {}).get("negative_control_pattern")
    if pattern:
        out["is_negative_control"] = out["library_id"].astype(str).str.contains(pattern, regex=True, na=False)
    else:
        out["is_negative_control"] = False

    out = out.drop_duplicates(subset=["library_id"]).reset_index(drop=True)
    return out[cols]


def _build_viruses(df: pd.DataFrame | None, filters: Dict[str, Any]) -> pd.DataFrame:
    cols = ["virus_id", "virus_family", "virus_name", "classification_evidence",
            "representative_accession", "mammal_associated", "eligible",
            "eligibility_reason"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)
    rules = filters.get("lineage_eligibility", {}) or {}
    allowed = {str(f).strip().lower() for f in rules.get("mammal_associated_families", [])}
    excluded = {str(f).strip().lower() for f in rules.get("exclude_families", [])}
    drop_unclassified = bool(rules.get("exclude_unclassified", True))

    out = pd.DataFrame(index=df.index)
    out["virus_id"] = df.get("virus_id").astype(str)
    family = df.get("virus_family")
    out["virus_family"] = family.fillna(UNKNOWN) if family is not None else UNKNOWN
    out["virus_name"] = df.get("virus_name")
    out["classification_evidence"] = df.get("classification_evidence")
    out["representative_accession"] = df.get("representative_accession")
    out["mammal_associated"] = df.get("mammal_associated")

    eligible, reason = [], []
    for fam in out["virus_family"].astype(str):
        key = fam.strip().lower()
        if key in excluded:
            eligible.append(False)
            reason.append("family excluded (non-mammalian host group)")
        elif key in {UNKNOWN, "", "nan", "unclassified"}:
            eligible.append(not drop_unclassified)
            reason.append("unclassified family")
        elif allowed and key not in allowed:
            eligible.append(False)
            reason.append("family not in the mammal-associated list")
        else:
            eligible.append(True)
            reason.append("eligible")
    out["eligible"] = eligible
    out["eligibility_reason"] = reason
    out = out.drop_duplicates(subset=["virus_id"]).reset_index(drop=True)
    return out[cols]


def _build_ecology(df: pd.DataFrame | None) -> pd.DataFrame:
    cols = ["site_id", "forest_cover", "humidity", "temperature", "precipitation", "diet"]
    if df is None or df.empty:
        return pd.DataFrame(columns=cols)
    out = pd.DataFrame(index=df.index)
    out["site_id"] = df.get("site_id").astype(str)
    for c in ("forest_cover", "humidity", "temperature", "precipitation"):
        out[c] = pd.to_numeric(df.get(c), errors="coerce") if c in df else np.nan
    out["diet"] = df.get("diet")
    return out.drop_duplicates(subset=["site_id"]).reset_index(drop=True)[cols]


# ---------------------------------------------------------------------------
# detection grid
# ---------------------------------------------------------------------------
def _build_detections(
    observed: pd.DataFrame | None,
    libraries: pd.DataFrame,
    viruses: pd.DataFrame,
    filters: Dict[str, Any],
    not_screened: pd.DataFrame | None = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Expand observed rows into a full library x lineage grid with states.

    The grid is the whole point: a model needs the zeros, but only the zeros
    that are real.
    """
    cols = ["library_id", "virus_id", "detection", "y", "mapped_reads",
            "coverage_breadth", "abundance", "evidence", "flag"]
    diagnostics: Dict[str, Any] = {}
    if libraries.empty or viruses.empty:
        return pd.DataFrame(columns=cols), {"reason": "no libraries or no viruses"}

    call = filters.get("detection_call", {}) or {}
    positive = call.get("positive", {}) or {}
    negative = call.get("negative_requires", {}) or {}
    min_reads = positive.get("min_mapped_reads")
    min_breadth = positive.get("min_coverage_breadth")
    min_abundance = positive.get("min_abundance")
    min_depth = negative.get("min_effective_depth")
    require_panel = bool(negative.get("screened_against_common_panel", True))

    lib_ids = libraries["library_id"].astype(str).tolist()
    vir_ids = viruses["virus_id"].astype(str).tolist()
    grid = pd.MultiIndex.from_product([lib_ids, vir_ids],
                                      names=["library_id", "virus_id"]).to_frame(index=False)

    obs = pd.DataFrame(columns=["library_id", "virus_id", "mapped_reads",
                                "coverage_breadth", "abundance"])
    if observed is not None and not observed.empty:
        obs = observed.copy()
        obs["library_id"] = obs["library_id"].astype(str)
        obs["virus_id"] = obs["virus_id"].astype(str)
        for c in ("mapped_reads", "coverage_breadth", "abundance"):
            obs[c] = pd.to_numeric(obs[c], errors="coerce") if c in obs else np.nan
        obs = obs.groupby(["library_id", "virus_id"], as_index=False).agg(
            mapped_reads=("mapped_reads", "max"),
            coverage_breadth=("coverage_breadth", "max"),
            abundance=("abundance", "max"),
        )

    df = grid.merge(obs, on=["library_id", "virus_id"], how="left")
    df["observed_row"] = df["mapped_reads"].notna() | df["coverage_breadth"].notna() | df["abundance"].notna()

    # Library-level screening status.
    panel = libraries.set_index("library_id")["screened_panel"].astype(str)
    depth = libraries.set_index("library_id")["effective_depth"]
    common_panel = panel.mode().iat[0] if len(panel) and panel.notna().any() else UNKNOWN
    lib_screened = panel.eq(common_panel) if require_panel else pd.Series(True, index=panel.index)
    if min_depth is not None:
        lib_screened = lib_screened & (depth.fillna(0) >= float(min_depth))

    df["library_screened"] = df["library_id"].map(lib_screened).fillna(False)

    # Cells explicitly recorded as not screened in the source deposit.
    if not_screened is not None and not not_screened.empty:
        ns = not_screened.copy()
        ns.columns = [c.lower() for c in ns.columns]
        lib_col = next((c for c in ns.columns if "lib" in c or c in {"sample", "run"}), ns.columns[0])
        vir_col = next((c for c in ns.columns if c not in {lib_col}), ns.columns[-1])
        ns_keys = set(zip(ns[lib_col].astype(str), ns[vir_col].astype(str)))
        df["explicitly_not_screened"] = [
            (lib, vir) in ns_keys for lib, vir in zip(df["library_id"], df["virus_id"])
        ]
    else:
        df["explicitly_not_screened"] = False

    # Positive call.
    is_positive = pd.Series(True, index=df.index)
    if min_reads is not None:
        is_positive &= df["mapped_reads"].fillna(-np.inf) >= float(min_reads)
    if min_breadth is not None:
        is_positive &= df["coverage_breadth"].fillna(-np.inf) >= float(min_breadth)
    if min_abundance is not None:
        is_positive &= df["abundance"].fillna(-np.inf) >= float(min_abundance)
    is_positive &= df["observed_row"]

    # A cell is assessable only if its library was screened against the common
    # panel at adequate depth, and the cell itself was not recorded as
    # unscreened. Crucially this applies to positives too: counting a positive
    # from an unscreened library while discarding that library's zeros would
    # bias detection upward exactly where sequencing effort was lowest.
    assessable = (
        df["library_screened"].to_numpy(dtype=bool)
        & (~df["explicitly_not_screened"].to_numpy(dtype=bool))
    )
    positive_mask = is_positive.to_numpy() & assessable

    detection = np.full(len(df), "unassessable", dtype=object)
    detection[positive_mask] = "positive"
    detection[(~is_positive.to_numpy()) & assessable] = "negative"
    df["detection"] = detection
    df["y"] = np.where(df["detection"] == "positive", 1.0,
                       np.where(df["detection"] == "negative", 0.0, np.nan))

    evidence = np.where(
        df["observed_row"], "published_table",
        np.where(df["detection"] == "negative", "screened_absent", "not_screened"),
    )
    # Keep the information that a detection exists outside the assessed
    # universe, so it is excluded from the likelihood but not lost.
    outside = is_positive.to_numpy() & (~assessable)
    evidence = np.where(outside, "positive_outside_assessed_universe", evidence)
    df["evidence"] = evidence
    df["flag"] = np.where(outside, "detected_but_library_not_screened", "")

    diagnostics = {
        "common_reference_panel": common_panel,
        "n_cells": int(len(df)),
        "n_positive": int((df["detection"] == "positive").sum()),
        "n_negative": int((df["detection"] == "negative").sum()),
        "n_unassessable": int((df["detection"] == "unassessable").sum()),
        "n_libraries_not_screened": int((~lib_screened).sum()),
        "n_explicit_not_screened_cells": int(df["explicitly_not_screened"].sum()),
        "observed_rows_supplied": int(len(obs)),
        "detection_rule": {
            "min_mapped_reads": min_reads,
            "min_coverage_breadth": min_breadth,
            "min_abundance": min_abundance,
            "negative_min_effective_depth": min_depth,
            "negative_requires_common_panel": require_panel,
        },
    }
    return df[cols], diagnostics


def _apply_contamination_filters(detections: pd.DataFrame, libraries: pd.DataFrame,
                                 filters: Dict[str, Any]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    rules = filters.get("contamination", {}) or {}
    diag: Dict[str, Any] = {}
    if detections.empty:
        return detections, diag
    df = detections.copy()
    batch = libraries.set_index("library_id").get("sequencing_batch")

    if rules.get("drop_negative_control_lineages") and "is_negative_control" in libraries:
        controls = set(
            libraries.loc[libraries["is_negative_control"].fillna(False), "library_id"].astype(str)
        )
        if controls:
            tainted = set(
                df.loc[df["library_id"].isin(controls) & (df["detection"] == "positive"), "virus_id"]
            )
            mask = df["virus_id"].isin(tainted) & (df["detection"] == "positive")
            df.loc[mask, "flag"] = (df.loc[mask, "flag"].astype(str) + ";detected_in_negative_control").str.strip(";")
            df.loc[mask, "detection"] = "unassessable"
            df.loc[mask, "y"] = np.nan
            diag["lineages_in_negative_controls"] = sorted(tainted)

    ratio = rules.get("index_hopping_ratio")
    if ratio and batch is not None and "abundance" in df:
        df["_batch"] = df["library_id"].map(batch)
        peak = df[df["detection"] == "positive"].groupby(["_batch", "virus_id"])["abundance"].transform("max")
        aligned = pd.Series(np.nan, index=df.index)
        aligned.loc[peak.index] = peak
        hop = (df["detection"] == "positive") & df["abundance"].notna() & aligned.notna() & (
            df["abundance"] < float(ratio) * aligned
        )
        df.loc[hop, "flag"] = (df.loc[hop, "flag"] + ";possible_index_hopping").str.strip(";")
        df.loc[hop, "detection"] = "negative"
        df.loc[hop, "y"] = 0.0
        diag["index_hopping_calls_downgraded"] = int(hop.sum())
        df = df.drop(columns=["_batch"])

    if rules.get("flag_single_batch_lineages") and batch is not None:
        pos = df[df["detection"] == "positive"].copy()
        pos["_batch"] = pos["library_id"].map(batch)
        per_lineage = pos.groupby("virus_id")["_batch"].nunique()
        single = set(per_lineage[per_lineage <= 1].index)
        mask = df["virus_id"].isin(single) & (df["detection"] == "positive")
        df.loc[mask, "flag"] = (df.loc[mask, "flag"] + ";single_batch_lineage").str.strip(";")
        diag["single_batch_lineages"] = sorted(single)

    return df, diag


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def build_canonical(raw_dir: str | Path, mapping: ResolvedMapping,
                    filters: Dict[str, Any]) -> Dict[str, Any]:
    """Build every canonical table plus a harmonization diagnostics record."""
    root = Path(raw_dir)
    loaded: Dict[str, pd.DataFrame | None] = {}
    crosswalk_parts: List[pd.DataFrame] = []
    for entity in ("sites", "hosts", "libraries", "viruses", "detections", "ecology"):
        df, _ = _load_entity(mapping, entity, root)
        loaded[entity] = df

    sites = _build_sites(loaded["sites"], mapping)
    hosts = _build_hosts(loaded["hosts"], mapping)
    libraries = _build_libraries(loaded["libraries"], mapping, filters)
    viruses = _build_viruses(loaded["viruses"], filters)
    ecology = _build_ecology(loaded["ecology"])

    # Libraries may carry their site only through their host pool.
    if not libraries.empty and not hosts.empty:
        missing_site = libraries["site_id"].isna() | (libraries["site_id"].astype(str) == "None")
        if missing_site.any():
            lookup = hosts.set_index("host_id")["site_id"]
            libraries.loc[missing_site, "site_id"] = libraries.loc[missing_site, "host_id"].map(lookup)

    not_screened = None
    ns_name = (mapping.extras or {}).get("not_screened_table")
    if ns_name:
        for candidate in root.rglob(Path(ns_name).name):
            not_screened = read_any_table(candidate)
            break

    detections, det_diag = _build_detections(
        loaded["detections"], libraries, viruses, filters, not_screened
    )
    detections, contam_diag = _apply_contamination_filters(detections, libraries, filters)

    for entity, frame, key in (
        ("sites", sites, "site_id"),
        ("hosts", hosts, "host_id"),
        ("libraries", libraries, "library_id"),
        ("viruses", viruses, "virus_id"),
    ):
        if frame.empty:
            continue
        res = mapping.entities.get(entity)
        source_col = None
        if res is not None:
            source_col = next((f.source_column for f in res.fields if f.canonical == key), None)
        crosswalk_parts.append(
            _crosswalk_rows(entity, frame[key], res.source_table if res else None,
                            source_col, "direct identifier carried through")
        )
    if not libraries.empty and "sra_run" in libraries:
        runs = libraries[["sra_run", "library_id"]].dropna()
        if len(runs):
            crosswalk_parts.append(
                pd.DataFrame(
                    {
                        "entity": "library_sra_run",
                        "source_id": runs["sra_run"].astype(str),
                        "canonical_id": runs["library_id"].astype(str),
                        "source_table": (mapping.entities["libraries"].source_table
                                         if "libraries" in mapping.entities else None),
                        "source_column": "sra_run",
                        "evidence": "library metadata links the SRA run to the library",
                        "provenance": "evonexus.harmonize.canonical",
                    }
                )
            )

    crosswalk = (
        pd.concat(crosswalk_parts, ignore_index=True)
        if crosswalk_parts
        else pd.DataFrame(columns=["entity", "source_id", "canonical_id", "source_table",
                                   "source_column", "evidence", "provenance"])
    ).drop_duplicates(subset=["entity", "source_id", "canonical_id"]).reset_index(drop=True)

    diagnostics = {
        "detections": det_diag,
        "contamination": contam_diag,
        "counts": {
            "sites": int(len(sites)),
            "hosts": int(len(hosts)),
            "libraries": int(len(libraries)),
            "viruses": int(len(viruses)),
            "viruses_eligible": int(viruses["eligible"].sum()) if "eligible" in viruses and len(viruses) else 0,
            "detections_cells": int(len(detections)),
            "ecology_sites": int(len(ecology)),
            "crosswalk_rows": int(len(crosswalk)),
        },
        "habitat_distribution": (
            sites["habitat"].value_counts().to_dict() if "habitat" in sites and len(sites) else {}
        ),
        "pool_status_distribution": (
            hosts["pool_status"].value_counts().to_dict() if "pool_status" in hosts and len(hosts) else {}
        ),
        "mixed_species_pools": int(hosts["mixed_species_pool"].fillna(False).sum()) if len(hosts) else 0,
        "effective_depth_definition": (
            libraries["effective_depth_definition"].iat[0] if len(libraries) else None
        ),
    }

    return {
        "tables": {
            "sites": sites,
            "hosts": hosts,
            "libraries": libraries,
            "viruses": viruses,
            "detections": detections,
            "ecology": ecology,
            "crosswalk": crosswalk,
        },
        "diagnostics": diagnostics,
    }
