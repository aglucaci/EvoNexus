// Data acquisition: manifests, downloads, inventory.
//
// Publishing contract, shared by every module: a process writes into a
// stage-named directory (`acquisition/`, `processed/`, `qc/`, …) and
// publishes to `params.outdir` itself. Naming the stage in both the output
// path and publishDir is what produced `results/acquisition/acquisition/`.

process ACQUIRE_DATA {
    label 'network'
    tag "${params.acquisition_mode}${params.synthetic ? ':synthetic' : ''}"
    publishDir "${params.outdir}", mode: params.publish_mode

    output:
    path 'acquisition/manifest.json',           emit: manifest
    path 'acquisition/inventory.json',          emit: inventory
    path 'acquisition/inventory.md',            emit: inventory_md
    path 'acquisition/raw',                     emit: raw
    path 'acquisition/run_record_acquire.json', emit: run_record
    // Only produced in synthetic mode; it carries the generating parameters
    // so the Bayesian stage can run a parameter-recovery check.
    path 'acquisition/raw/synthetic/synthetic_truth.json', emit: truth, optional: true

    script:
    def synthetic = params.synthetic ? '--synthetic' : ''
    def offline   = params.offline   ? '--offline'   : ''
    """
    python -m evonexus.cli acquire \\
        --out-dir acquisition \\
        --mode ${params.acquisition_mode} \\
        --seed ${params.seed} \\
        ${synthetic} ${offline}
    """
}

process DOWNLOAD_READS {
    label 'network'
    tag "${run_accession}"
    storeDir "${params.datadir}/raw/sra/${run_accession}"

    input:
    val run_accession

    output:
    tuple val(run_accession), path("${run_accession}*.fastq.gz"), emit: reads
    path "${run_accession}.md5",                                  emit: checksums

    when:
    params.acquisition_mode == 'raw_reads'

    script:
    """
    prefetch --max-size u ${run_accession}
    fasterq-dump --split-files --threads ${task.cpus} --outdir . ${run_accession}
    pigz -p ${task.cpus} ${run_accession}*.fastq || gzip ${run_accession}*.fastq
    md5sum ${run_accession}*.fastq.gz > ${run_accession}.md5
    """
}
