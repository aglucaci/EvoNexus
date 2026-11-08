# Changelog

All notable changes to EvoNexus are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — 2026-11-08

First working pipeline.

### Added
- Nextflow DSL2 pipeline with `test`, `laptop`, `processed`, `full`,
  `evolution` and `manuscript` profiles, plus `docker`, `singularity`, `conda`
  and `hpc` execution profiles.
- Acquisition of the BtCN-Virome Figshare deposit, SRA RunInfo for
  PRJNA862267, PRJNA997461 and PRJNA929070, and GenBank homologs, with a
  permanent manifest (identifiers, URLs, retrieval dates, checksums, versions).
- Source inventory reporting real table headers, dtypes and missingness, and
  stating which pipeline prerequisites the deposit actually meets.
- Metadata harmonization to a canonical schema with a crosswalk, refusing to
  run until a human approves the column mapping.
- Seven mandatory data gates (G1–G7); a critical failure stops the primary
  habitat analysis and emits an explanatory QC report.
- Library × lineage response matrix distinguishing positives, assay-tested
  negatives and unassessable cells.
- Hierarchical Bayesian detection model M0–M5 with PyMC/NUTS, plus a
  self-contained Pólya-Gamma Gibbs engine with interweaved centered and
  non-centered scale updates for environments without PyMC.
- Posterior standardization of the settlement-versus-wild estimand,
  per-family and per-host strata, and sensitivity across priors, lineage
  resolutions and sample restrictions.
- Posterior cross-host connectivity with a marginal null (no host × family
  interaction) and a sampling-preserving permutation null.
- Held-out site and region cross-validation that integrates over unseen
  random effects, with calibration curves and proper scoring rules.
- Evolutionary subworkflow (MAFFT, IQ-TREE 2, HyPhy) with eligibility gates
  that produce documented exclusions rather than misleading numbers.
- Automated Figures 1–5, Tables 1–4 and Supplementary 1–5 with generated
  legends; unavailable analyses are emitted as "NOT GENERATED".
- 59 unit and integration tests, and GitHub Actions CI.
