// Manuscript production: figures, source tables, legends and the run report.

process GENERATE_REPORT {
    label 'light'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path results

    output:
    path 'figures/*',    emit: figures,  optional: true
    path 'tables/*',     emit: tables,   optional: true
    path 'manuscript/*', emit: manuscript

    when:
    !params.skip_report

    script:
    """
    python -m evonexus.cli report --results-dir ${results} --out-dir .
    """
}

process RENDER_DATA_DICTIONARY {
    label 'tiny'
    publishDir "${params.outdir}/docs", mode: params.publish_mode

    output:
    path 'DATA_DICTIONARY.md', emit: dictionary

    script:
    """
    python -m evonexus.cli schema --markdown DATA_DICTIONARY.md
    """
}
