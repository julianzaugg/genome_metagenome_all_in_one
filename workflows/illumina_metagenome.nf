/*
 * ILLUMINA_METAGENOME — placeholder; full implementation lands in Task 4.
 */

include { INPUT_CHECK } from '../subworkflows/local/input_check'

workflow ILLUMINA_METAGENOME {
    if (!params.input) { error "Mode 'illumina_metagenome' requires --input <samplesheet.csv>" }
    INPUT_CHECK(params.input, params.mode)
    log.warn "[gmaio] illumina_metagenome placeholder — being implemented."
}
