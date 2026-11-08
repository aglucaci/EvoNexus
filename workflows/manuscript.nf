include { GENERATE_REPORT; RENDER_DATA_DICTIONARY } from '../modules/reporting/main.nf'

workflow MANUSCRIPT {
    take:
    results

    main:
    RENDER_DATA_DICTIONARY()
    GENERATE_REPORT(results)

    emit:
    manuscript = GENERATE_REPORT.out.manuscript
}
