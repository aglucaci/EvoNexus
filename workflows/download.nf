include { ACQUIRE_DATA; DOWNLOAD_READS } from '../modules/acquisition/main.nf'

workflow DOWNLOAD {
    main:
    // ACQUIRE_DATA produces all three acquisition outputs itself: the
    // permanent manifest, the unmodified raw files, and the inventory of what
    // the deposit actually contains.
    ACQUIRE_DATA()

    // Raw reads only when the acquisition mode asks for them. Run accessions
    // come from the RunInfo tables the acquisition step already recorded.
    if (params.acquisition_mode == 'raw_reads') {
        def runs = ACQUIRE_DATA.out.manifest
            .splitJson(path: 'sources')
            .filter { it instanceof Map && it.fastq_download }
            .flatMap { it.fastq_download.run_accessions ?: [] }
            .unique()
        DOWNLOAD_READS(runs)
    }

    emit:
    manifest  = ACQUIRE_DATA.out.manifest
    inventory = ACQUIRE_DATA.out.inventory
    raw       = ACQUIRE_DATA.out.raw
    truth     = ACQUIRE_DATA.out.truth
}
