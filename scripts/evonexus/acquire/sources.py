"""Clients for the public repositories EvoNexus reads from.

Only the standard library is used for HTTP so that acquisition works in a
minimal container. Every request is retried with exponential backoff and every
downloaded file is checksummed by the caller.
"""

from __future__ import annotations

import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import pandas as pd

from evonexus.core.io import sha256_file

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SRA_RUNINFO = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi"


class AcquisitionError(RuntimeError):
    """Raised when a source cannot be retrieved after retries."""


@dataclass
class HttpClient:
    user_agent: str = "EvoNexus/0.1 (reproducible research pipeline)"
    max_retries: int = 5
    backoff_seconds: float = 2.0
    timeout_seconds: float = 120.0
    offline: bool = False

    def _request(self, url: str, data: bytes | None = None) -> bytes:
        if self.offline:
            raise AcquisitionError(f"offline mode: refusing to fetch {url}")
        last: Exception | None = None
        for attempt in range(self.max_retries):
            try:
                req = urllib.request.Request(
                    url, data=data, headers={"User-Agent": self.user_agent}
                )
                with urllib.request.urlopen(req, timeout=self.timeout_seconds) as resp:
                    return resp.read()
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
                last = exc
                if attempt == self.max_retries - 1:
                    break
                time.sleep(self.backoff_seconds * (2 ** attempt))
        raise AcquisitionError(f"failed to fetch {url}: {last}")

    def get_bytes(self, url: str) -> bytes:
        return self._request(url)

    def get_json(self, url: str) -> Any:
        return json.loads(self._request(url).decode("utf-8"))

    def get_text(self, url: str) -> str:
        return self._request(url).decode("utf-8", errors="replace")

    def download(self, url: str, dest: str | Path) -> Dict[str, Any]:
        """Download ``url`` to ``dest`` and return a manifest entry."""
        target = Path(dest)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.get_bytes(url)
        target.write_bytes(payload)
        return {
            "url": url,
            "path": str(target),
            "bytes": target.stat().st_size,
            "sha256": sha256_file(target),
        }


# ---------------------------------------------------------------------------
# Figshare
# ---------------------------------------------------------------------------
def figshare_file_list(client: HttpClient, article_id: int | str,
                       api_endpoint: str | None = None) -> List[Dict[str, Any]]:
    """List the files in a Figshare article without downloading them."""
    endpoint = api_endpoint or f"https://api.figshare.com/v2/articles/{article_id}"
    meta = client.get_json(endpoint)
    files = meta.get("files", [])
    out = []
    for f in files:
        out.append(
            {
                "name": f.get("name"),
                "download_url": f.get("download_url"),
                "size_bytes": f.get("size"),
                "supplied_md5": f.get("computed_md5") or f.get("supplied_md5"),
                "figshare_file_id": f.get("id"),
            }
        )
    return out


def figshare_article_metadata(client: HttpClient, article_id: int | str,
                              api_endpoint: str | None = None) -> Dict[str, Any]:
    endpoint = api_endpoint or f"https://api.figshare.com/v2/articles/{article_id}"
    meta = client.get_json(endpoint)
    return {
        "title": meta.get("title"),
        "doi": meta.get("doi"),
        "version": meta.get("version"),
        "published_date": meta.get("published_date"),
        "license": (meta.get("license") or {}).get("name"),
        "n_files": len(meta.get("files", [])),
    }


# ---------------------------------------------------------------------------
# NCBI SRA
# ---------------------------------------------------------------------------
def sra_runinfo(client: HttpClient, bioproject: str) -> pd.DataFrame:
    """Fetch the SRA RunInfo table for a BioProject accession.

    Returns an empty dataframe (rather than raising) when the project has no
    runs, so that an absent secondary dataset degrades to a recorded gap.
    """
    search = (
        f"{EUTILS}/esearch.fcgi?db=sra&term={urllib.parse.quote(bioproject)}"
        "[BioProject]&retmax=10000&usehistory=y&retmode=json"
    )
    result = client.get_json(search).get("esearchresult", {})
    webenv = result.get("webenv")
    query_key = result.get("querykey")
    count = int(result.get("count", 0))
    if not count or not webenv:
        return pd.DataFrame()
    fetch = (
        f"{EUTILS}/efetch.fcgi?db=sra&rettype=runinfo&retmode=text"
        f"&WebEnv={webenv}&query_key={query_key}&retmax={count}"
    )
    text = client.get_text(fetch)
    if not text.strip():
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO(text))
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
    return df


# ---------------------------------------------------------------------------
# GenBank / nuccore
# ---------------------------------------------------------------------------
def genbank_fetch(client: HttpClient, accessions: Iterable[str],
                  batch_size: int = 200, rettype: str = "fasta") -> str:
    """Fetch sequence records for a list of accessions, batched."""
    accs = [a.strip() for a in accessions if a and a.strip()]
    chunks: List[str] = []
    for i in range(0, len(accs), batch_size):
        batch = accs[i : i + batch_size]
        url = (
            f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(batch)}"
            f"&rettype={rettype}&retmode=text"
        )
        chunks.append(client.get_text(url))
        time.sleep(0.4)  # stay inside NCBI's unauthenticated rate limit
    return "\n".join(chunks)


def genbank_search(client: HttpClient, query: str, retmax: int = 500) -> List[str]:
    url = (
        f"{EUTILS}/esearch.fcgi?db=nuccore&term={urllib.parse.quote(query)}"
        f"&retmax={retmax}&retmode=json"
    )
    result = client.get_json(url).get("esearchresult", {})
    return list(result.get("idlist", []))
