// Data acquisition: manifests, downloads, inventory.

process ACQUIRE_DATA {
    label 'network'
    tag "${params.acquisition_mode}${params.synthetic ? ':synthetic' : ''}"
    publishDir "${params.outdir}/acquisition", mode: params.publish_mode

    output:
    path 'acquisition/manifest.json',   emit: manifest
    path 'acquisition/inventory.json',  emit: inventory
    path 'acquisition/inventory.md',    emit: inventory_md
    path 'acquisition/raw',             emit: raw
    path 'acquisition/run_record_acquire.json', emit: run_record

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

process INVENTORY_SOURCES {
    label 'light'
    publishDir "${params.outdir}/qc", mode: params.publish_mode

    input:
    path raw

    output:
    path 'inventory/inventory.json', emit: inventory
    path 'inventory/inventory.md',   emit: report

    script:
    """
    python -m evonexus.cli inventory --raw-dir ${raw} --out-dir inventory
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
