// Bayesian inference, connectivity and predictive validation.

process FIT_MODEL {
    label 'heavy'
    tag "${params.prior_variant}:${params.engine}"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(resolution), path(frame)
    path processed
    path host_covariance
    path truth
    val  primary_allowed

    output:
    path 'bayesian', emit: bayesian

    script:
    def models    = params.models.join(' ')
    def synthetic = params.synthetic ? '--synthetic' : ''
    def sens      = params.run_sensitivity ? '' : '--no-sensitivity'
    // A critical QC gate failure routes to the descriptive branch rather than
    // aborting: the CLI then writes unadjusted rates and states plainly that
    // no habitat effect was estimated.
    def skip      = (primary_allowed.toString().toLowerCase() == 'true') ? '' : '--skip-primary'
    def hostcov   = host_covariance.name != 'NO_FILE' ? "--host-covariance ${host_covariance}" : ''
    def truth_arg = truth.name != 'NO_FILE' ? "--truth ${truth}" : ''
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
        ${hostcov} ${truth_arg} ${synthetic} ${sens} ${skip}
    """
}

process COMPUTE_CONNECTIVITY {
    label 'standard'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(resolution), path(frame)
    path bayesian
    path host_covariance

    output:
    path 'connectivity', emit: connectivity

    script:
    def synthetic = params.synthetic ? '--synthetic' : ''
    def hostcov   = host_covariance.name != 'NO_FILE' ? "--host-covariance ${host_covariance}" : ''
    """
    python -m evonexus.cli connectivity \\
        --posterior-dir ${bayesian} \\
        --model ${params.primary_model} \\
        --frame ${frame} \\
        --out-dir connectivity \\
        --permutations ${params.permutations} \\
        --seed ${params.seed} ${hostcov} ${synthetic}
    """
}

process CROSS_VALIDATE {
    label 'heavy'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(resolution), path(frame)
    path bayesian

    output:
    path 'validation', emit: validation

    script:
    def groups    = params.cv_group_columns.join(' ')
    def synthetic = params.synthetic ? '--synthetic' : ''
    """
    python -m evonexus.cli validate \\
        --frame ${frame} \\
        --out-dir validation \\
        --posterior-dir ${bayesian} \\
        --model ${params.primary_model} \\
        --group-columns ${groups} \\
        --folds ${params.cv_folds} \\
        --engine ${params.engine} \\
        --draws ${params.draws} --tune ${params.tune} --chains ${params.chains} \\
        ${synthetic}
    """
}
