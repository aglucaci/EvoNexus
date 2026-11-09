include { CROSS_VALIDATE } from '../modules/statistics/main.nf'

workflow VALIDATION {
    take:
    frame
    bayesian

    main:
    CROSS_VALIDATE(frame, bayesian)

    emit:
    validation = CROSS_VALIDATE.out.validation
}
