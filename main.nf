#!/usr/bin/env nextflow
/*
 * EvoNexus — evolutionary and ecological networks of viral host sharing.
 *
 * The pipeline is deliberately gated. Data provenance and sample
 * harmonization are validated before any Bayesian model runs, because an
 * incorrect relationship between an SRA run, a sequencing library, a bat
 * species and a viral lineage would invalidate every downstream estimate.
 *
 *   discovery -> harmonize -> QC gates -> [pass] -> inference -> manuscript
 *                                      \-> [fail] -> descriptive only + report
 */

nextflow.enable.dsl = 2

// Includes must be at script level: Nextflow's parser rejects an `include`
// inside a workflow body. The host-phylogeny process is therefore always
// imported and conditionally *invoked* below.
include { BUILD_HOST_PHYLOGENY } from './modules/phylogenetics/main.nf'
include { DOWNLOAD }     from './workflows/download.nf'
include { PREPROCESS }   from './workflows/preprocess.nf'
include { BAYESIAN }     from './workflows/bayesian.nf'
include { CONNECTIVITY } from './workflows/connectivity.nf'
include { VALIDATION }   from './workflows/validation.nf'
include { EVOLUTION }    from './workflows/evolution.nf'
include { MANUSCRIPT }   from './workflows/manuscript.nf'

def helpMessage() {
    log.info """
    ===========================================================================
    EvoNexus ${workflow.manifest.version}
    ===========================================================================

    Usage:
      nextflow run . -profile test
      nextflow run . -profile laptop
      nextflow run . -profile processed --acquisition_mode processed_only
      nextflow run . -profile full,singularity --acquisition_mode raw_reads
      nextflow run . -profile evolution --sequences data/external/viral.fasta
      nextflow run . -profile manuscript

    Profiles:
      test        Deterministic synthetic data; validates the whole DAG.
      laptop      Small public subset; limited chains and viral panel.
      processed   Full published processed dataset (default analysis).
      full        All eligible public reads and metadata; HPC via SLURM.
      evolution   Molecular evolutionary analyses on selected lineages.
      manuscript  Regenerate figures, tables and legends from existing results.
      docker / singularity / conda / hpc may be combined with the above.

    Key parameters:
      --acquisition_mode        processed_only (default) | raw_reads
      --require_mapping_approval  true (default). Harmonization refuses to run
                                  until config/metadata_mapping.yaml is approved.
      --fail_on_critical_gate   false (default). On a critical QC failure the
                                primary habitat model is skipped and a
                                descriptive, clearly-labelled branch runs.
      --models                  Model ladder members to fit (default M0..M4)
      --engine                  auto | pymc_nuts | polya_gamma_gibbs
      --draws --tune --chains --seed
      --sequences               Viral FASTA for the evolution workflow
      --host_newick             Bat phylogeny enabling model M5
      --outdir                  Results directory

    Outputs map to the manuscript as:
      Figure 1  results/figures/Figure1.*   sampling structure
      Figure 2  results/figures/Figure2.*   habitat inference
      Figure 3  results/figures/Figure3.*   cross-host connectivity
      Figure 4  results/figures/Figure4.*   viral evolutionary structure
      Figure 5  results/figures/Figure5.*   validation
    ===========================================================================
    """.stripIndent()
}

workflow {

    // `return` is not permitted in a workflow body under Nextflow's strict
    // syntax, so the three operating modes are expressed as one if/else chain
    // rather than as early exits.
    if (params.help) {
        helpMessage()
    }
    else if (params.mode == 'manuscript') {
        // Regenerate figures, tables and legends from existing results.
        // No acquisition, no sampling.
        log.info "EvoNexus ${workflow.manifest.version} · manuscript mode · outdir: ${params.outdir}"
        MANUSCRIPT(Channel.fromPath(params.outdir, type: 'dir', checkIfExists: true))
    }
    else {
        log.info """
        EvoNexus ${workflow.manifest.version}  ·  mode=${params.mode}  ·  \
synthetic=${params.synthetic}  ·  engine=${params.engine}
        outdir: ${params.outdir}
        """.stripIndent()

        if (params.synthetic) {
            log.warn 'SYNTHETIC MODE — outputs validate software only and are ' +
                     'watermarked. They are not scientific results.'
        }

        def no_file = file("${projectDir}/assets/NO_FILE")
        def sequences_ch   = params.sequences   ? Channel.fromPath(params.sequences,   checkIfExists: true)
                                                : Channel.value(no_file)
        def host_newick_ch = params.host_newick ? Channel.fromPath(params.host_newick, checkIfExists: true)
                                                : Channel.value(no_file)

        // ---- 01 acquire and inventory ------------------------------------
        DOWNLOAD()

        // ---- 02-05 harmonize, curate, build the matrix, run the gates -----
        PREPROCESS(DOWNLOAD.out.inventory, DOWNLOAD.out.raw)

        // A host phylogeny is optional. Without it model M5 is skipped and the
        // connectivity-versus-relatedness panel is marked as not generated.
        def host_cov_ch = Channel.value(no_file)
        if (params.host_newick) {
            BUILD_HOST_PHYLOGENY(PREPROCESS.out.processed, host_newick_ch)
            host_cov_ch = BUILD_HOST_PHYLOGENY.out.covariance.ifEmpty(no_file)
        }

        // Synthetic runs carry their generating parameters, which enables the
        // parameter-recovery check. Real runs have none, and it is skipped.
        def truth_ch = DOWNLOAD.out.truth.ifEmpty(no_file)

        // ---- 06 Bayesian inference ----------------------------------------
        BAYESIAN(
            PREPROCESS.out.primary_frame,
            PREPROCESS.out.processed,
            host_cov_ch,
            truth_ch,
            PREPROCESS.out.primary_allowed
        )

        // ---- 07-10 connectivity, validation, evolution ---------------------
        CONNECTIVITY(PREPROCESS.out.primary_frame, BAYESIAN.out.bayesian, host_cov_ch)
        VALIDATION(PREPROCESS.out.primary_frame, BAYESIAN.out.bayesian)

        // Connectivity results feed the evolutionary integration step. When
        // the primary model was skipped there are none, and EVOLUTION runs
        // without them rather than waiting forever.
        def connectivity_table = CONNECTIVITY.out.connectivity
            .map { dir -> dir.list().find { it.startsWith('connectivity.') }
                              ? file("${dir}/" + dir.list().find { it.startsWith('connectivity.') })
                              : no_file }
            .ifEmpty(no_file)

        EVOLUTION(PREPROCESS.out.processed, sequences_ch, connectivity_table, host_newick_ch)

        // ---- 11-12 manuscript ----------------------------------------------
        // Reporting reads the published results tree, so it waits for every
        // analysis to finish and is then handed that directory. Returning a
        // Channel from inside `map` would nest channels; emit the path itself.
        def results_ch = PREPROCESS.out.qc
            .mix(BAYESIAN.out.bayesian,
                 CONNECTIVITY.out.connectivity.ifEmpty(no_file),
                 VALIDATION.out.validation.ifEmpty(no_file),
                 EVOLUTION.out.evolution)
            .collect()
            .map { file(params.outdir) }

        MANUSCRIPT(results_ch)
    }
}
