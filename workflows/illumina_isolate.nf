/*
 * ILLUMINA_ISOLATE — SCAFFOLD (v1: wiring + TODO stubs)
 *
 * Target flow (per-sample then cross-sample within each `group`):
 *   fastp QC -> shovill assembly
 *     -> CheckM1/2 + GTDB-Tk (genome_taxonomy_qc)
 *     -> bakta + MLST + AMRFinderPlus + ISEScan (isolate_annotation)
 *     -> geNomad -> CheckV -> ANI cluster (mobile_elements)
 *     -> CoverM mapping
 *   group-wise (collect()/groupTuple on meta.group):
 *     -> fastANI all-vs-all (comparative)
 *     -> chewBBACA cgMLST
 *     -> parsnp -> gubbins (phylogenomics)
 *     -> panaroo -> core-gene tree (pangenome)
 */

include { INPUT_CHECK } from '../subworkflows/local/input_check'

workflow ILLUMINA_ISOLATE {

    if (!params.input) { error "Mode 'illumina_isolate' requires --input <samplesheet.csv>" }

    INPUT_CHECK(params.input, params.mode)

    // Channels: INPUT_CHECK.out.reads_short / .host / .meta
    // TODO: SHORT_READ_QC(fastp)
    // TODO: ASSEMBLY_ILLUMINA_ISO(shovill)
    // TODO: GENOME_TAXONOMY_QC(gtdbtk + checkm1/2)
    // TODO: ISOLATE_ANNOTATION(bakta + mlst + amrfinderplus + isescan)
    // TODO: MOBILE_ELEMENTS(genomad -> checkv -> ANI cluster)
    // TODO: READ_MAPPING(coverm)
    // group-wise:
    // TODO: COMPARATIVE(fastani)  ; CGMLST(chewbacca)
    // TODO: PHYLOGENOMICS(parsnp -> gubbins)  ; PANGENOME(panaroo -> tree)

    log.warn "[gmaio] mode 'illumina_isolate' is a SCAFFOLD — channel wiring is in place but processes are not yet implemented."
}
