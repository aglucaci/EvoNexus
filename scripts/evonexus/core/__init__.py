"""Shared infrastructure: configuration, tabular IO, schema, provenance."""

from evonexus.core.config import load_config, load_yaml, resolve_variant
from evonexus.core.io import (
    TABLE_SUFFIX,
    read_table,
    write_table,
    write_json,
    read_json,
    sha256_file,
    sha256_bytes,
)
from evonexus.core.provenance import RunRecord, record_environment
from evonexus.core.schema import CANONICAL_SCHEMAS, validate_table, SchemaError

__all__ = [
    "load_config",
    "load_yaml",
    "resolve_variant",
    "TABLE_SUFFIX",
    "read_table",
    "write_table",
    "write_json",
    "read_json",
    "sha256_file",
    "sha256_bytes",
    "RunRecord",
    "record_environment",
    "CANONICAL_SCHEMAS",
    "validate_table",
    "SchemaError",
]
