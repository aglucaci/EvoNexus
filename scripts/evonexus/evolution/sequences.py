"""FASTA handling and lineage selection for the evolutionary subworkflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd


def read_fasta(path: str | Path) -> Dict[str, str]:
    """Minimal FASTA reader; keeps the full header line as the key."""
    records: Dict[str, str] = {}
    name: str | None = None
    chunks: List[str] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].strip()
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def write_fasta(records: Dict[str, str], path: str | Path, line_width: int = 60) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), line_width):
                handle.write(seq[i: i + line_width] + "\n")
    return target


def ungapped_length(seq: str) -> int:
    return sum(1 for c in seq if c not in "-.~ ")


def mean_pairwise_divergence(records: Dict[str, str], max_pairs: int = 5000,
                             seed: int = 0) -> float:
    """Uncorrected p-distance averaged over aligned pairs."""
    seqs = [s.upper() for s in records.values()]
    n = len(seqs)
    if n < 2:
        return float("nan")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        return float("nan")      # not an alignment
    rng = np.random.default_rng(seed)
    pairs: List[Tuple[int, int]] = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if len(pairs) > max_pairs:
        idx = rng.choice(len(pairs), max_pairs, replace=False)
        pairs = [pairs[i] for i in idx]
    distances = []
    for i, j in pairs:
        a, b = seqs[i], seqs[j]
        comparable = [(x, y) for x, y in zip(a, b) if x in "ACGTU" and y in "ACGTU"]
        if len(comparable) < 50:
            continue
        diff = sum(1 for x, y in comparable if x != y)
        distances.append(diff / len(comparable))
    return float(np.mean(distances)) if distances else float("nan")


def gap_fractions(records: Dict[str, str]) -> Tuple[Dict[str, float], np.ndarray]:
    seqs = list(records.values())
    if not seqs:
        return {}, np.zeros(0)
    length = len(seqs[0])
    per_seq = {
        name: sum(1 for c in seq if c in "-.~") / max(len(seq), 1)
        for name, seq in records.items()
    }
    matrix = np.array([[1 if c in "-.~" else 0 for c in seq] for seq in seqs])
    per_column = matrix.mean(axis=0) if matrix.size else np.zeros(length)
    return per_seq, per_column


def select_lineages(
    viruses: pd.DataFrame,
    detections: pd.DataFrame,
    sequences: Dict[str, str],
    config: Dict[str, Any],
) -> Dict[str, Any]:
    """Pick lineages with enough homologous sequence for phylogenetics.

    Lineages that do not qualify are returned in ``excluded`` with the reason,
    so the manuscript can state what was left out and why.
    """
    rules = config.get("sequence_selection", {}) or {}
    min_seqs = int(rules.get("min_sequences_per_lineage", 8))
    min_len = int(rules.get("min_sequence_length", 300))
    max_seqs = int(rules.get("max_sequences_per_lineage", 500))
    focal = set(rules.get("focal_families") or [])

    # Sequences are assigned to lineages by a prefix in the FASTA header.
    by_lineage: Dict[str, Dict[str, str]] = {}
    for header, seq in sequences.items():
        token = header.split()[0] if header.split() else header
        lineage = token.split("|")[0]
        by_lineage.setdefault(lineage, {})[header] = seq

    selected: Dict[str, Dict[str, str]] = {}
    excluded: List[Dict[str, Any]] = []
    rows: List[Dict[str, Any]] = []

    families = (
        viruses.set_index("virus_id")["virus_family"].astype(str).to_dict()
        if len(viruses) else {}
    )
    detected_hosts: Dict[str, int] = {}
    if len(detections):
        pos = detections[detections["detection"] == "positive"]
        detected_hosts = pos.groupby("virus_id")["library_id"].nunique().to_dict()

    for lineage, records in by_lineage.items():
        family = families.get(lineage, "unknown")
        kept = {h: s for h, s in records.items() if ungapped_length(s) >= min_len}
        reason = None
        if focal and family not in focal:
            reason = f"family {family} is not in the focal list"
        elif len(kept) < min_seqs:
            reason = f"{len(kept)} sequences >= {min_len} nt (minimum {min_seqs})"
        if reason:
            excluded.append({"lineage": lineage, "family": family, "reason": reason,
                             "n_sequences": len(records)})
            continue
        if len(kept) > max_seqs:
            kept = dict(list(kept.items())[:max_seqs])
        selected[lineage] = kept
        rows.append(
            {
                "lineage": lineage,
                "family": family,
                "n_sequences": len(kept),
                "median_length": int(np.median([ungapped_length(s) for s in kept.values()])),
                "n_libraries_detected": int(detected_hosts.get(lineage, 0)),
            }
        )

    return {
        "selected": selected,
        "table": pd.DataFrame(rows),
        "excluded": pd.DataFrame(excluded),
        "n_selected": len(selected),
        "n_excluded": len(excluded),
    }
