include { ACQUIRE_DATA; INVENTORY_SOURCES; DOWNLOAD_READS } from '../modules/acquisition/main.nf'

workflow DOWNLOAD {
    main:
    ACQUIRE_DATA()
    INVENTORY_SOURCES(ACQUIRE_DATA.out.raw)

    // Raw reads only when the acquisition mode asks for them. Run accessions
    // come from the RunInfo tables the acquisition step already recorded.
    reads = Channel.empty()
    if (params.acquisition_mode == 'raw_reads') {
        runs = ACQUIRE_DATA.out.manifest
            .splitJson()
            .flatMap { entry -> entry.value instanceof List ? entry.value : [] }
            .filter { it instanceof Map && it.fastq_download }
            .flatMap { it.fastq_download.run_accessions ?: [] }
            .unique()
        reads = DOWNLOAD_READS(runs).reads
    }

    emit:
    manifest  = ACQUIRE_DATA.out.manifest
    inventory = INVENTORY_SOURCES.out.inventory
    raw       = ACQUIRE_DATA.out.raw
    reads     = reads
}
