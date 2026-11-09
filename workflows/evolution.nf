include { RUN_EVOLUTION } from '../modules/phylogenetics/main.nf'

workflow EVOLUTION {
    take:
    processed
    sequences
    connectivity
    host_newick

    main:
    RUN_EVOLUTION(processed, sequences, connectivity, host_newick)

    emit:
    evolution = RUN_EVOLUTION.out.evolution
}
