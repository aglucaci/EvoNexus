include { CROSS_VALIDATE } from '../modules/statistics/main.nf'

workflow VALIDATION {
    take:
    frame
    posteriors

    main:
    CROSS_VALIDATE(frame, posteriors)

    emit:
    all      = CROSS_VALIDATE.out.all.collect()
    metadata = CROSS_VALIDATE.out.metadata
}
