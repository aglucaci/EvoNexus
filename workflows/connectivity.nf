include { COMPUTE_CONNECTIVITY } from '../modules/statistics/main.nf'

workflow CONNECTIVITY {
    take:
    frame
    posteriors
    host_covariance

    main:
    COMPUTE_CONNECTIVITY(frame, posteriors, host_covariance)

    emit:
    all      = COMPUTE_CONNECTIVITY.out.all.collect()
    metadata = COMPUTE_CONNECTIVITY.out.metadata
}
