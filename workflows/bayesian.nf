include { FIT_MODEL } from '../modules/statistics/main.nf'

workflow BAYESIAN {
    take:
    frame
    processed
    host_covariance
    truth
    primary_allowed

    main:
    FIT_MODEL(frame, processed, host_covariance, truth, primary_allowed)

    emit:
    bayesian = FIT_MODEL.out.bayesian
}
