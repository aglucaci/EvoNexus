// Viral evolutionary analyses. Each stage may end in a documented exclusion.

process BUILD_HOST_PHYLOGENY {
    label 'light'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path processed
    path newick

    output:
    path 'host_phylogeny/host_covariance.json', emit: covariance, optional: true

    script:
    """
    python -m evonexus.cli evolution \\
        --processed-dir ${processed} \\
        --out-dir host_phylogeny \\
        --host-newick ${newick}
    """
}

process RUN_EVOLUTION {
    label 'heavy'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path processed
    path sequences
    path connectivity
    path newick

    output:
    path 'evolution', emit: evolution

    script:
    def seqs      = sequences.name    != 'NO_FILE' ? "--sequences ${sequences}" : ''
    def tree      = newick.name       != 'NO_FILE' ? "--host-newick ${newick}"  : ''
    def conn      = connectivity.name != 'NO_FILE' ? "--connectivity ${connectivity}" : ''
    def synthetic = params.synthetic ? '--synthetic' : ''
    """
    python -m evonexus.cli evolution \\
        --processed-dir ${processed} \\
        --out-dir evolution \\
        ${seqs} ${tree} ${conn} ${synthetic}
    """
}
