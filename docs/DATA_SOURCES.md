# Data sources

Every source is declared in [`config/datasets.yaml`](../config/datasets.yaml).
Nothing is downloaded implicitly: the acquisition step resolves each source to
a concrete file list, records a manifest with checksums, and only then
retrieves it.

## Primary dataset — BtCN-Virome

**Yan et al. (2025).** *A compendium of 8,176 bat RNA viral metagenomes
reveals ecological drivers and circulation dynamics.*

The study reports 4,143 bats, 40 bat species, 52 locations and **363
sequencing libraries**. The 363 libraries — not the number of animals and not
the number of assembled viral sequences — define the observation structure for
library-level modelling. Verify these figures against the publication before
submission; the pipeline reports whatever the files actually contain.

| Data | Source | Retrieval |
| --- | --- | --- |
| Processed viral catalog | Figshare [10.6084/m9.figshare.24515056](https://doi.org/10.6084/m9.figshare.24515056) | Figshare API file manifest, then per-file download |
| Ecological metadata | Same deposit | Discovered by inventory, not assumed |
| Published alignments and trees | Same deposit | Retrieved and audited |
| Raw metatranscriptomes | SRA [PRJNA862267](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA862267), [PRJNA997461](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA997461) | Entrez RunInfo; FASTQ only in `raw_reads` mode |
| Viral reference sequences | GenBank accessions from the publication | NCBI E-utilities, batched |

**The pipeline does not assume the deposit contains a complete library ×
lineage count matrix.** The inventory inspects what is there and reports
whether a matrix can be built from the processed files or whether raw-read
mapping is required (`expectations.verdict` in `inventory.json`).

## Independent dataset — individual-bat sequencing

| Data | Source | Role |
| --- | --- | --- |
| Individual-bat metatranscriptomes | SRA [PRJNA929070](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA929070) | Secondary evaluation |

This dataset is used to ask whether host-sharing summaries remain
interpretable without pooled sequencing. **It is not an independent
replication of the settlement-versus-wild habitat contrast**, and the pipeline
never presents it as one.

## Environmental covariates

Opt-in and disabled by default. ERA5-Land and land-cover products are joined
only when the source metadata carries geographic and temporal resolution that
supports them; place pre-extracted per-site covariates under
`data/external/`. Covariates are included in model M4 only if at least 80% of
modelled cells have a value.

## Acquisition modes

- `processed_only` (default) — processed tables, metadata, alignments, trees,
  plus SRA RunInfo (small, and it is what links runs to libraries).
- `raw_reads` — additionally retrieves FASTQ via `prefetch`/`fasterq-dump`
  with per-run checksums and resumable caching, for reconstructing
  library-level detections against a consistent reference catalog.

## Licensing

The pipeline code is MIT. **The datasets are not.** Check the license recorded
in `manifest.json` for the Figshare deposit and the terms attached to each
NCBI record before redistributing any derived data. Cite the source
publications independently of citing EvoNexus.

## What the manifest records

For every source: repository identifier, filenames, URLs, retrieval date in
UTC, SHA-256 checksums, deposit version and license, and — for sources that
could not be reached — the failure, so a partial run is visible rather than
silent.
