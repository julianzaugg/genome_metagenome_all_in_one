/*
 * NANOPORE_ISOLATE — SCAFFOLD (v1: wiring + TODO stubs)
 *
 * Target flow (per-sample then cross-sample within each `group`):
 *   pod5 -> dorado basecall/demux (optional) -> porechop -> fastplong (long_read_qc)
 *     -> autocycler assembly -> dorado polish
 *        -> polypolish  [ONLY if meta.has_short_reads — hybrid]
 *        -> dnaapler reorientation
 *     -> CheckM1/2 + GTDB-Tk (genome_taxonomy_qc)
 *     -> bakta + MLST + AMRFinderPlus + ISEScan (isolate_annotation)
 *     -> geNomad -> CheckV -> ANI cluster (mobile_elements)
 *     -> CoverM mapping (long + short where available)
 *   group-wise: fastANI ; chewBBACA cgMLST ; parsnp -> gubbins ; panaroo -> tree
 */

include { INPUT_CHECK } from '../subworkflows/local/input_check'

workflow NANOPORE_ISOLATE {

    if (!params.input) { error "Mode 'nanopore_isolate' requires --input <samplesheet.csv>" }

    INPUT_CHECK(params.input, params.mode)

    // Hybrid note: INPUT_CHECK.out.reads_short carries only samples that ALSO have
    // Illumina reads; join it onto the polished assemblies to gate POLYPOLISH per-sample.
    //
    // TODO: (optional) DORADO_BASECALLER + DORADO_DEMUX from INPUT_CHECK.out.pod5
    // TODO: LONG_READ_QC(porechop -> fastplong)
    // TODO: ASSEMBLY_NANOPORE_ISO(autocycler -> dorado polish -> [polypolish if short] -> dnaapler)
    // TODO: GENOME_TAXONOMY_QC ; ISOLATE_ANNOTATION ; MOBILE_ELEMENTS ; READ_MAPPING
    // group-wise: COMPARATIVE ; CGMLST ; PHYLOGENOMICS ; PANGENOME

    log.warn "[gmaio] mode 'nanopore_isolate' is a SCAFFOLD — channel wiring is in place but processes are not yet implemented."
}
