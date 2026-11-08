"""Driver for the ACQUIRE_DATA step.

Three outputs, always:

* ``manifest.json`` — permanent record of repository ids, filenames, URLs,
  retrieval dates, checksums and versions.
* a raw directory holding source files **unmodified**.
* ``inventory.json`` / ``inventory.md`` — what the deposit actually contains.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from evonexus.acquire import sources
from evonexus.acquire.inventory import inventory_directory, write_inventory
from evonexus.core.config import load_config
from evonexus.core.io import sha256_file, write_json
from evonexus.core.provenance import RunRecord, utc_now

MODES = ("processed_only", "raw_reads")


def _client(cfg: Dict[str, Any], offline: bool) -> sources.HttpClient:
    http = cfg.get("http", {}) or {}
    return sources.HttpClient(
        user_agent=http.get("user_agent", "EvoNexus/0.1"),
        max_retries=int(http.get("max_retries", 5)),
        backoff_seconds=float(http.get("backoff_seconds", 2.0)),
        timeout_seconds=float(http.get("timeout_seconds", 120)),
        offline=offline,
    )


def _acquire_figshare(client, name, spec, raw_dir, record) -> Dict[str, Any]:
    entry: Dict[str, Any] = {"source": name, "kind": "figshare", "status": "pending"}
    dest = Path(raw_dir) / name
    try:
        meta = sources.figshare_article_metadata(
            client, spec.get("article_id"), spec.get("api_endpoint")
        )
        files = sources.figshare_file_list(
            client, spec.get("article_id"), spec.get("api_endpoint")
        )
    except sources.AcquisitionError as exc:
        record.warn(f"{name}: could not reach Figshare ({exc})")
        entry.update(status="unavailable", error=str(exc))
        return entry

    entry["article_metadata"] = meta
    entry["files_listed"] = files
    downloaded: List[Dict[str, Any]] = []
    for f in files:
        url = f.get("download_url")
        if not url:
            continue
        target = dest / (f.get("name") or f"file_{f.get('figshare_file_id')}")
        if target.exists():
            downloaded.append(
                {"url": url, "path": str(target), "bytes": target.stat().st_size,
                 "sha256": sha256_file(target), "cached": True}
            )
            continue
        try:
            got = client.download(url, target)
            got["cached"] = False
            got["supplied_md5"] = f.get("supplied_md5")
            downloaded.append(got)
            record.note(f"{name}: downloaded {target.name} ({got['bytes'] / 1e6:.1f} MB)")
        except sources.AcquisitionError as exc:
            record.warn(f"{name}: failed to download {url}: {exc}")
            downloaded.append({"url": url, "error": str(exc)})
    entry["downloaded"] = downloaded
    entry["status"] = "ok" if any("sha256" in d for d in downloaded) else "empty"
    return entry


def _acquire_sra(client, name, spec, raw_dir, record, mode) -> Dict[str, Any]:
    accession = spec.get("accession")
    entry: Dict[str, Any] = {"source": name, "kind": "sra_bioproject",
                             "accession": accession, "status": "pending"}
    dest = Path(raw_dir) / name
    dest.mkdir(parents=True, exist_ok=True)
    runinfo_path = dest / f"{accession}_runinfo.csv"
    if runinfo_path.exists():
        df = pd.read_csv(runinfo_path)
        entry["cached"] = True
    else:
        try:
            df = sources.sra_runinfo(client, accession)
        except sources.AcquisitionError as exc:
            record.warn(f"{name}: RunInfo unavailable ({exc})")
            entry.update(status="unavailable", error=str(exc))
            return entry
        df.to_csv(runinfo_path, index=False)
    entry["runinfo_path"] = str(runinfo_path)
    entry["runinfo_sha256"] = sha256_file(runinfo_path) if runinfo_path.exists() else None
    entry["n_runs"] = int(len(df))
    if len(df) and "bases" in df.columns:
        entry["total_bases"] = int(pd.to_numeric(df["bases"], errors="coerce").fillna(0).sum())
    record.note(f"{name}: {len(df)} SRA runs recorded for {accession}")

    if mode == "raw_reads":
        entry["fastq_download"] = {
            "performed": False,
            "reason": (
                "FASTQ retrieval is delegated to the Nextflow DOWNLOAD_READS "
                "process, which uses prefetch/fasterq-dump with per-run "
                "checksum verification and resumable caching."
            ),
            "run_accessions": df["Run"].astype(str).tolist() if "Run" in df.columns else [],
        }
    entry["status"] = "ok"
    return entry


def _acquire_genbank(client, name, spec, raw_dir, record) -> Dict[str, Any]:
    entry: Dict[str, Any] = {"source": name, "kind": "ncbi_nuccore", "status": "pending"}
    query_file = Path(spec.get("query_file", "config/genbank_queries.txt"))
    if not query_file.exists():
        entry.update(status="skipped", reason=f"no query file at {query_file}")
        record.note(f"{name}: skipped — no accession/query list at {query_file}")
        return entry
    queries = [
        line.strip()
        for line in query_file.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not queries:
        entry.update(status="skipped", reason="query file empty")
        return entry
    dest = Path(raw_dir) / name
    dest.mkdir(parents=True, exist_ok=True)
    out_path = dest / "homologs.fasta"
    try:
        fasta = sources.genbank_fetch(client, queries, int(spec.get("batch_size", 200)))
    except sources.AcquisitionError as exc:
        record.warn(f"{name}: GenBank fetch failed ({exc})")
        entry.update(status="unavailable", error=str(exc))
        return entry
    out_path.write_text(fasta, encoding="utf-8")
    entry.update(
        status="ok",
        path=str(out_path),
        sha256=sha256_file(out_path),
        n_queries=len(queries),
        n_records=fasta.count(">"),
    )
    return entry


def run_acquire(
    out_dir: str | Path,
    mode: str = "processed_only",
    datasets_config: str = "datasets",
    offline: bool = False,
    synthetic: bool = False,
    synthetic_seed: int = 20250101,
    synthetic_libraries: int = 120,
    synthetic_viruses: int = 30,
) -> Dict[str, Any]:
    """Acquire sources (or generate synthetic ones) and inventory the result."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")

    out = Path(out_dir)
    raw_dir = out / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    record = RunRecord(step="acquire", params={"mode": mode, "offline": offline,
                                               "synthetic": synthetic})
    record.synthetic = synthetic

    entries: List[Dict[str, Any]] = []

    if synthetic:
        from evonexus.simulate import simulate

        record.note("synthetic mode: generating deterministic test data")
        sim = simulate(
            raw_dir / "synthetic",
            n_libraries=synthetic_libraries,
            n_viruses=synthetic_viruses,
            seed=synthetic_seed,
        )
        entries.append(
            {
                "source": "synthetic",
                "kind": "simulation",
                "status": "ok",
                "seed": synthetic_seed,
                "files": sim["paths"],
                "warning": "SYNTHETIC DATA — software validation only, not a scientific result.",
            }
        )
    else:
        cfg = load_config(datasets_config)
        client = _client(cfg, offline)
        for name, spec in (cfg.get("sources") or {}).items():
            if not spec.get("enabled", True):
                entries.append({"source": name, "status": "disabled"})
                continue
            tier = spec.get("tier")
            if tier == "reads" and mode == "processed_only":
                # Still record RunInfo: it is small and it is what links SRA
                # runs to libraries, which the harmonization step needs.
                entries.append(_acquire_sra(client, name, spec, raw_dir, record, "processed_only"))
                continue
            kind = spec.get("kind")
            if kind == "figshare":
                entries.append(_acquire_figshare(client, name, spec, raw_dir, record))
            elif kind == "sra_bioproject":
                entries.append(_acquire_sra(client, name, spec, raw_dir, record, mode))
            elif kind == "ncbi_nuccore":
                entries.append(_acquire_genbank(client, name, spec, raw_dir, record))
            else:
                entries.append({"source": name, "status": "unsupported_kind", "kind": kind})

    manifest = {
        "retrieved_at_utc": utc_now(),
        "mode": mode,
        "synthetic": synthetic,
        "offline": offline,
        "sources": entries,
    }
    manifest_path = write_json(manifest, out / "manifest.json")
    record.add_output(manifest_path, role="manifest")

    inv = inventory_directory(raw_dir)
    inv["synthetic"] = synthetic
    inv_paths = write_inventory(inv, out)
    for p in inv_paths.values():
        record.add_output(p, role="inventory")

    failures = [e for e in entries if e.get("status") in {"unavailable", "empty"}]
    status = "ok" if not failures else "partial"
    if failures:
        record.warn(
            f"{len(failures)} source(s) unavailable: "
            + ", ".join(str(e.get("source")) for e in failures)
        )
    record.finish(out / "run_record_acquire.json", status=status)

    return {
        "manifest": str(manifest_path),
        "inventory_json": str(inv_paths["json"]),
        "inventory_md": str(inv_paths["markdown"]),
        "raw_dir": str(raw_dir),
        "status": status,
        "n_sources": len(entries),
    }
