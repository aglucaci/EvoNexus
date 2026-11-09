include { HARMONIZE_METADATA; RUN_QC_GATES } from '../modules/metadata/main.nf'
include { CREATE_MATRIX }                    from '../modules/viromics/main.nf'

workflow PREPROCESS {
    take:
    inventory
    raw

    main:
    // Harmonization refuses to run until the column mapping is approved and
    // every required field resolves, so this is the pipeline's real gate.
    HARMONIZE_METADATA(inventory, raw)

    RUN_QC_GATES(HARMONIZE_METADATA.out.processed)

    CREATE_MATRIX(HARMONIZE_METADATA.out.processed, params.resolutions as List)

    emit:
    processed       = HARMONIZE_METADATA.out.processed
    qc              = RUN_QC_GATES.out.qc
    // The gate verdict routes the run: true fits the primary habitat model,
    // false falls through to the descriptive branch. Read from a file rather
    // than an `env` output, which older Nextflow releases do not support.
    primary_allowed = RUN_QC_GATES.out.qc.map { qc ->
                          file("${qc}/primary_model_allowed.txt").text.trim() }
    frames          = CREATE_MATRIX.out.frame
    // The first configured resolution drives modelling; the rest are
    // sensitivity analyses. Emitted as an expression rather than a local:
    // a `def` inside the workflow body is not visible to the emit block.
    primary_frame   = CREATE_MATRIX.out.frame.filter { resolution, _frame ->
                          resolution == params.resolutions[0] }
}
