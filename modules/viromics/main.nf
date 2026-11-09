// The library x lineage response matrix.
//
// There is no separate reference-building process: the curated viral catalog,
// with its eligibility flags and exclusion reasons, is produced by
// harmonization as `processed/viruses.*`. A process that only copied it under
// another name would add a failure mode without adding information.

process CREATE_MATRIX {
    label 'light'
    tag "${resolution}"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path processed
    each resolution

    output:
    tuple val(resolution), path("matrix/analysis_frame_${resolution}.*"), emit: frame
    path "matrix/*",                                                      emit: all

    script:
    def synthetic = params.synthetic ? '--synthetic' : ''
    """
    python -m evonexus.cli matrix \\
        --processed-dir ${processed} \\
        --out-dir matrix \\
        --resolution ${resolution} ${synthetic}
    """
}
