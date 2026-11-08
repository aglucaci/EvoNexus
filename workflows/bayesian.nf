include { FIT_MODEL } from '../modules/statistics/main.nf'

workflow BAYESIAN {
    take:
    frame
    processed
    host_covariance
    primary_allowed

    main:
    FIT_MODEL(frame, processed, host_covariance, primary_allowed)

    emit:
    all        = FIT_MODEL.out.all.collect()
    status     = FIT_MODEL.out.status
    posteriors = FIT_MODEL.out.all.collect()
}
