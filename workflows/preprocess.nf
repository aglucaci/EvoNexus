include { HARMONIZE_METADATA; RUN_QC_GATES } from '../modules/metadata/main.nf'
include { BUILD_REFERENCE; CREATE_MATRIX }  from '../modules/viromics/main.nf'

workflow PREPROCESS {
    take:
    inventory
    raw

    main:
    // Harmonization refuses to run until the column mapping is approved and
    // every required field resolves, so this is the pipeline's real gate.
    HARMONIZE_METADATA(inventory, raw)
    processed = HARMONIZE_METADATA.out.all.collect()

    RUN_QC_GATES(processed)
    BUILD_REFERENCE(processed)

    resolutions = Channel.fromList(params.resolutions as List)
    CREATE_MATRIX(processed, resolutions)

    // The primary resolution drives modelling; the others are sensitivity.
    primary_frame = CREATE_MATRIX.out.frame.filter { res, _f -> res == params.resolutions[0] }

    emit:
    processed       = processed
    qc_gates        = RUN_QC_GATES.out.gates
    qc_report       = RUN_QC_GATES.out.report
    primary_allowed = RUN_QC_GATES.out.primary_allowed
    frames          = CREATE_MATRIX.out.frame
    primary_frame   = primary_frame
}
