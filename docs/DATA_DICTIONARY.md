# Data dictionary

Generated from `evonexus.core.schema`. Do not edit by hand.

## `sites.parquet`

- Unit of observation: sampling site
- Primary key: `site_id`

| Column | Type | Required | Nullable | Allowed values | Description |
| --- | --- | --- | --- | --- | --- |
| `site_id` | string | yes | no | — | Canonical site identifier |
| `region` | string | no | yes | — | — |
| `country` | string | no | yes | — | — |
| `habitat` | category | yes | yes | `settlement`, `wild`, `unknown` | Settlement / wild / unknown; unknown is excluded from the contrast |
| `habitat_raw` | string | no | yes | — | Verbatim source value |
| `latitude` | float | no | yes | — | — |
| `longitude` | float | no | yes | — | — |
| `geo_resolution` | string | no | yes | — | Stated precision of the coordinates |

## `hosts.parquet`

- Unit of observation: animal or documented pool
- Primary key: `host_id`
- Foreign keys: `site_id` → `sites.site_id`

| Column | Type | Required | Nullable | Allowed values | Description |
| --- | --- | --- | --- | --- | --- |
| `host_id` | string | yes | no | — | — |
| `site_id` | string | yes | yes | — | — |
| `host_species` | string | yes | yes | — | Verbatim species, or 'mixed' / 'unknown' |
| `host_species_resolved` | string | yes | yes | — | Harmonized species name used for modelling |
| `host_genus` | string | no | yes | — | — |
| `host_family` | string | no | yes | — | — |
| `pool_status` | category | yes | yes | `individual`, `pool`, `unknown` | — |
| `pool_size` | float | no | yes | — | Number of animals; NaN when undocumented |
| `mixed_species_pool` | bool | no | yes | — | True when the pool spans more than one species |
| `collection_date` | string | no | yes | — | — |

## `libraries.parquet`

- Unit of observation: sequencing library
- Primary key: `library_id`
- Foreign keys: `host_id` → `hosts.host_id`, `site_id` → `sites.site_id`

| Column | Type | Required | Nullable | Allowed values | Description |
| --- | --- | --- | --- | --- | --- |
| `library_id` | string | yes | no | — | — |
| `host_id` | string | yes | yes | — | — |
| `site_id` | string | yes | yes | — | — |
| `tissue` | string | yes | yes | — | — |
| `tissue_group` | string | yes | yes | — | Controlled tissue grouping |
| `read_depth` | float | yes | yes | — | Total reads reported by the source |
| `effective_depth` | float | yes | yes | — | Depth denominator used in the model; see viral_filters.yaml |
| `effective_depth_definition` | string | yes | yes | — | — |
| `sequencing_batch` | string | no | yes | — | — |
| `sra_run` | string | no | yes | — | — |
| `bioproject` | string | no | yes | — | — |
| `collection_date` | string | no | yes | — | — |
| `collection_time_numeric` | float | no | yes | — | Decimal years since the first collection date |
| `screened_panel` | string | no | yes | — | Identifier of the reference panel this library was screened against |
| `is_negative_control` | bool | no | yes | — | — |

## `viruses.parquet`

- Unit of observation: viral lineage
- Primary key: `virus_id`

| Column | Type | Required | Nullable | Allowed values | Description |
| --- | --- | --- | --- | --- | --- |
| `virus_id` | string | yes | no | — | — |
| `virus_family` | string | yes | yes | — | — |
| `virus_name` | string | no | yes | — | — |
| `classification_evidence` | string | no | yes | — | — |
| `representative_accession` | string | no | yes | — | — |
| `mammal_associated` | bool | no | yes | — | — |
| `eligible` | bool | no | yes | — | Passes viral_filters.yaml lineage eligibility |
| `eligibility_reason` | string | no | yes | — | — |

## `detections.parquet`

- Unit of observation: library x viral lineage
- Primary key: `library_id`, `virus_id`
- Foreign keys: `library_id` → `libraries.library_id`, `virus_id` → `viruses.virus_id`

| Column | Type | Required | Nullable | Allowed values | Description |
| --- | --- | --- | --- | --- | --- |
| `library_id` | string | yes | no | — | — |
| `virus_id` | string | yes | no | — | — |
| `detection` | category | yes | no | `positive`, `negative`, `unassessable` | positive / negative / unassessable — never silently zero |
| `y` | float | no | yes | — | 1.0, 0.0 or NaN mirroring `detection` |
| `mapped_reads` | float | no | yes | — | — |
| `coverage_breadth` | float | no | yes | — | — |
| `abundance` | float | no | yes | — | — |
| `evidence` | string | no | yes | — | published_table | read_mapping | not_screened |
| `flag` | string | no | yes | — | Contamination / index-hopping flags, semicolon separated |

## `ecology.parquet`

- Unit of observation: site
- Primary key: `site_id`
- Foreign keys: `site_id` → `sites.site_id`

| Column | Type | Required | Nullable | Allowed values | Description |
| --- | --- | --- | --- | --- | --- |
| `site_id` | string | yes | no | — | — |
| `forest_cover` | float | no | yes | — | — |
| `humidity` | float | no | yes | — | — |
| `temperature` | float | no | yes | — | — |
| `precipitation` | float | no | yes | — | — |
| `diet` | string | no | yes | — | — |

## `crosswalk.parquet`

- Unit of observation: identifier relationship
- Primary key: `entity`, `source_id`, `canonical_id`

| Column | Type | Required | Nullable | Allowed values | Description |
| --- | --- | --- | --- | --- | --- |
| `entity` | string | yes | no | — | — |
| `source_id` | string | yes | no | — | — |
| `canonical_id` | string | yes | no | — | — |
| `source_table` | string | no | yes | — | — |
| `source_column` | string | no | yes | — | — |
| `evidence` | string | no | yes | — | — |
| `provenance` | string | no | yes | — | — |
