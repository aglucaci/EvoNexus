// Metadata harmonization and the mandatory data gates.

process HARMONIZE_METADATA {
    label 'light'
    publishDir "${params.outdir}/processed", mode: params.publish_mode, pattern: 'processed/*'
    publishDir "${params.outdir}/qc",        mode: params.publish_mode, pattern: 'processed/mapping_report.md'

    input:
    path inventory
    path raw

    output:
    path 'processed/*',                     emit: all
    path 'processed/libraries*',            emit: libraries
    path 'processed/detections*',           emit: detections
    path 'processed/harmonize_diagnostics.json', emit: diagnostics
    path 'processed/mapping_report.md',     emit: mapping_report

    script:
    def approval  = params.require_mapping_approval ? '' : '--skip-approval'
    def guess     = params.allow_role_guess ? '--allow-role-guess' : ''
    def synthetic = params.synthetic ? '--synthetic' : ''
    def overrides = params.mapping_overrides ? "--overrides ${params.mapping_overrides}" : ''
    """
    python -m evonexus.cli harmonize \\
        --inventory ${inventory} \\
        --raw-dir ${raw} \\
        --out-dir processed \\
        ${approval} ${guess} ${synthetic} ${overrides}
    """
}

process RUN_QC_GATES {
    label 'light'
    publishDir "${params.outdir}/qc", mode: params.publish_mode

    input:
    path processed

    output:
    path 'qc/qc_gates.json',           emit: gates
    path 'qc/qc_report.md',            emit: report
    path 'qc/exclusion_waterfall*',    emit: waterfall
    path 'qc/qc_gate_summary*',        emit: summary
    env  PRIMARY_ALLOWED,              emit: primary_allowed

    script:
    def synthetic = params.synthetic ? '--synthetic' : ''
    def strict    = params.fail_on_critical_gate ? '--fail-on-critical' : ''
    """
    python -m evonexus.cli qc \\
        --processed-dir ${processed} \\
        --out-dir qc ${synthetic} ${strict}

    PRIMARY_ALLOWED=\$(python -c "import json;print(json.load(open('qc/qc_gates.json'))['primary_model_allowed'])")
    """
}
