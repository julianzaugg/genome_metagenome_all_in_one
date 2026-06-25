/*
 * NANOPORE_METAGENOME — SCAFFOLD (v1: wiring + TODO stubs)
 *
 * Target flow:
 *   pod5 -> dorado basecall/demux (optional) -> porechop -> fastplong (long_read_qc)
 *     -> sylph + singlem (read_profiling) -> host_removal
 *     -> myloasm assembly -> dorado polish (assembly_nanopore_meta)
 *     -> Aviary binning -> CoverM dereplication -> CoverM mapping
 *     -> gene_catalogue (pyrodigal -> cd-hit -> DRAM)   [shared with illumina_metagenome]
 *     -> GTDB-Tk + CheckM1/2 + nonpareil + genomespot (genome_taxonomy_qc)
 *     -> geNomad -> CheckV -> ANI cluster (mobile_elements)
 */

include { INPUT_CHECK } from '../subworkflows/local/input_check'

workflow NANOPORE_METAGENOME {

    if (!params.input) { error "Mode 'nanopore_metagenome' requires --input <samplesheet.csv>" }

    INPUT_CHECK(params.input, params.mode)

    // Channels available: INPUT_CHECK.out.reads_long / .pod5 / .reads_short (hybrid) / .host / .meta
    // TODO: LONG_READ_QC(porechop -> fastplong)
    // TODO: READ_PROFILING(sylph + singlem)
    // TODO: HOST_REMOVAL
    // TODO: ASSEMBLY_NANOPORE_META(myloasm >=0.5.1 -> dorado polish)
    // TODO: BINNING(aviary) ; BIN_DEREPLICATION(coverm cluster) ; READ_MAPPING(coverm)
    // TODO: GENE_CATALOGUE(pyrodigal -> cdhit -> dram)   <-- reuse subworkflow from illumina_metagenome
    // TODO: GENOME_TAXONOMY_QC(gtdbtk + checkm1/2 + nonpareil + genomespot)
    // TODO: MOBILE_ELEMENTS(genomad -> checkv -> ANI cluster)

    log.warn "[gmaio] mode 'nanopore_metagenome' is a SCAFFOLD — channel wiring is in place but processes are not yet implemented."
}
