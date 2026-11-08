// Bayesian inference, connectivity and predictive validation.

process FIT_MODEL {
    label 'heavy'
    tag "${params.prior_variant}:${params.engine}"
    publishDir "${params.outdir}/bayesian", mode: params.publish_mode

    input:
    tuple val(resolution), path(frame)
    path processed
    path host_covariance
    val  primary_allowed

    output:
    path 'bayesian/*',                          emit: all
    path 'bayesian/bayesian_status.json',       emit: status
    path 'bayesian/habitat_contrast*',          emit: contrast,   optional: true
    path 'bayesian/model_comparison*',          emit: comparison, optional: true
    path 'bayesian/posterior_*',                emit: posteriors, optional: true

    script:
    def models    = params.models.join(' ')
    def synthetic = params.synthetic ? '--synthetic' : ''
    def sens      = params.run_sensitivity ? '' : '--no-sensitivity'
    def skip      = (primary_allowed.toString().toLowerCase() == 'true') ? '' : '--skip-primary'
    def hostcov   = host_covariance.name != 'NO_FILE' ? "--host-covariance ${host_covariance}" : ''
    def truth     = params.synthetic ? "--truth \$(ls ${params.outdir}/acquisition/raw/synthetic/synthetic_truth.json 2>/dev/null || echo '')" : ''
    """
    python -m evonexus.cli bayesian \\
        --frame ${frame} \\
        --out-dir bayesian \\
        --models ${models} \\
        --prior-variant ${params.prior_variant} \\
        --engine ${params.engine} \\
        --draws ${params.draws} --tune ${params.tune} --chains ${params.chains} \\
        --seed ${params.seed} \\
        --processed-dir ${processed} \\
        ${hostcov} ${synthetic} ${sens} ${skip} ${truth}
    """
}

process COMPUTE_CONNECTIVITY {
    label 'standard'
    publishDir "${params.outdir}/connectivity", mode: params.publish_mode

    input:
    tuple val(resolution), path(frame)
    path posterior_dir
    path host_covariance

    output:
    path 'connectivity/*', emit: all
    path 'connectivity/connectivity_metadata.json', emit: metadata

    script:
    def synthetic = params.synthetic ? '--synthetic' : ''
    def hostcov   = host_covariance.name != 'NO_FILE' ? "--host-covariance ${host_covariance}" : ''
    """
    python -m evonexus.cli connectivity \\
        --posterior-dir ${posterior_dir} \\
        --model ${params.primary_model} \\
        --frame ${frame} \\
        --out-dir connectivity \\
        --permutations ${params.permutations} \\
        --seed ${params.seed} ${hostcov} ${synthetic}
    """
}

process CROSS_VALIDATE {
    label 'heavy'
    publishDir "${params.outdir}/validation", mode: params.publish_mode

    input:
    tuple val(resolution), path(frame)
    path posterior_dir

    output:
    path 'validation/*', emit: all
    path 'validation/validation_metadata.json', emit: metadata

    script:
    def groups    = params.cv_group_columns.join(' ')
    def synthetic = params.synthetic ? '--synthetic' : ''
    """
    python -m evonexus.cli validate \\
        --frame ${frame} \\
        --out-dir validation \\
        --posterior-dir ${posterior_dir} \\
        --model ${params.primary_model} \\
        --group-columns ${groups} \\
        --folds ${params.cv_folds} \\
        --engine ${params.engine} \\
        --draws ${params.draws} --tune ${params.tune} --chains ${params.chains} \\
        ${synthetic}
    """
}
