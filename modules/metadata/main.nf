// Metadata harmonization and the mandatory data gates.

process HARMONIZE_METADATA {
    label 'light'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path inventory
    path raw

    output:
    // The whole directory is emitted, because every downstream stage takes a
    // `--processed-dir` and needs the canonical tables together rather than
    // as a loose list of files.
    path 'processed', emit: processed

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
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path processed

    output:
    path 'qc',                  emit: qc
    env('PRIMARY_ALLOWED'),     emit: primary_allowed

    script:
    def synthetic = params.synthetic ? '--synthetic' : ''
    def strict    = params.fail_on_critical_gate ? '--fail-on-critical' : ''
    """
    python -m evonexus.cli qc \\
        --processed-dir ${processed} \\
        --out-dir qc ${synthetic} ${strict}

    # Read back the verdict so the workflow can route to the primary model or
    # to the descriptive branch. A critical gate failure is not an error here:
    # it is a documented outcome.
    PRIMARY_ALLOWED=\$(python -c "import json;print(json.load(open('qc/qc_gates.json'))['primary_model_allowed'])")
    """
}
