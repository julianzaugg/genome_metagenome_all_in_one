#!/usr/bin/env nextflow
/*
 * gmaio — genome & metagenome all-in-one
 * Entry point. Routes on --mode to one of four track workflows (or DB download).
 */

nextflow.enable.dsl = 2

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

include { ILLUMINA_METAGENOME } from './workflows/illumina_metagenome'
include { NANOPORE_METAGENOME } from './workflows/nanopore_metagenome'
include { ILLUMINA_ISOLATE    } from './workflows/illumina_isolate'
include { NANOPORE_ISOLATE    } from './workflows/nanopore_isolate'
include { DOWNLOAD_DBS         } from './workflows/download_dbs'

def VALID_MODES = ['illumina_metagenome', 'nanopore_metagenome',
                   'illumina_isolate', 'nanopore_isolate', 'download_dbs']

workflow {

    // --- Help ---
    if (params.help) {
        log.info paramsHelp("nextflow run . -profile bunya --mode illumina_metagenome --input samplesheet.csv")
        exit 0
    }

    // --- Mode validation ---
    if (!params.mode) {
        error "You must supply --mode <${VALID_MODES.join(' | ')}>"
    }
    if (!VALID_MODES.contains(params.mode)) {
        error "Unknown --mode '${params.mode}'. Valid: ${VALID_MODES.join(', ')}"
    }

    // --- Param validation (nf-schema) ---
    if (params.validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    // --- Route ---
    switch (params.mode) {
        case 'download_dbs':        DOWNLOAD_DBS();         break
        case 'illumina_metagenome': ILLUMINA_METAGENOME();  break
        case 'nanopore_metagenome': NANOPORE_METAGENOME();  break
        case 'illumina_isolate':    ILLUMINA_ISOLATE();     break
        case 'nanopore_isolate':    NANOPORE_ISOLATE();     break
    }
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n[gmaio] Pipeline completed successfully (mode=${params.mode}). Results in: ${params.outdir}\n"
        : "\n[gmaio] Pipeline completed with errors (mode=${params.mode}).\n" )
}
