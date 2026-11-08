// Viral reference handling and the library x lineage response matrix.

process BUILD_REFERENCE {
    label 'light'
    publishDir "${params.outdir}/viromics", mode: params.publish_mode

    input:
    path processed

    output:
    path 'reference/viral_catalog*', emit: catalog, optional: true

    script:
    """
    mkdir -p reference
    cp -f ${processed}/viruses* reference/ 2>/dev/null || true
    for f in reference/viruses*; do [ -e "\$f" ] && mv "\$f" "\${f/viruses/viral_catalog}"; done
    """
}

process CREATE_MATRIX {
    label 'light'
    tag "${resolution}"
    publishDir "${params.outdir}/matrix", mode: params.publish_mode

    input:
    path processed
    val resolution

    output:
    tuple val(resolution), path("matrix/analysis_frame_${resolution}*"), emit: frame
    path "matrix/detection_matrix_${resolution}*",                        emit: matrix
    path "matrix/matrix_metadata_${resolution}.json",                     emit: metadata
    path "matrix/matrix_ledger_${resolution}*",                           emit: ledger

    script:
    def synthetic = params.synthetic ? '--synthetic' : ''
    """
    python -m evonexus.cli matrix \\
        --processed-dir ${processed} \\
        --out-dir matrix \\
        --resolution ${resolution} ${synthetic}
    """
}
