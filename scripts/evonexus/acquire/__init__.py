"""Data acquisition: manifests, downloads and source inventory."""

from evonexus.acquire.inventory import inventory_directory, write_inventory
from evonexus.acquire.run import run_acquire
from evonexus.acquire.sources import (
    HttpClient,
    figshare_file_list,
    genbank_fetch,
    sra_runinfo,
)

__all__ = [
    "HttpClient",
    "figshare_file_list",
    "sra_runinfo",
    "genbank_fetch",
    "inventory_directory",
    "write_inventory",
    "run_acquire",
]
