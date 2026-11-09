include { COMPUTE_CONNECTIVITY } from '../modules/statistics/main.nf'

workflow CONNECTIVITY {
    take:
    frame
    bayesian
    host_covariance

    main:
    COMPUTE_CONNECTIVITY(frame, bayesian, host_covariance)

    emit:
    connectivity = COMPUTE_CONNECTIVITY.out.connectivity
}
