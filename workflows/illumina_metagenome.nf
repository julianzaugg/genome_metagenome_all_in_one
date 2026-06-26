/*
 * ILLUMINA_METAGENOME — full end-to-end workflow.
 *
 *  fastp QC ─┬─> sylph + singlem (read profiling, on raw reads)
 *            └─> host removal ─> metaSPAdes ─┬─> pyrodigal ─> gene catalogue ─> DRAM
 *                                            ├─> geNomad ─> CheckV ─> ANI cluster
 *                                            └─> Aviary binning ─> CheckM2
 *                                                  ─> CoverM dereplication
 *                                                       ─> CoverM mapping (per sample)
 *                                                       ─> GTDB-Tk + CheckM1 + GenomeSPOT
 *  reads ─> nonpareil
 *
 * Cross-sample steps (catalogue, dereplication, GTDB-Tk, geNomad pooling) gather
 * per-sample outputs with .collect().
 */

include { INPUT_CHECK }        from '../subworkflows/local/input_check'
include { READ_PROFILING }     from '../subworkflows/local/read_profiling'
include { HOST_REMOVAL }       from '../subworkflows/local/host_removal'
include { GENE_CATALOGUE }     from '../subworkflows/local/gene_catalogue'
include { GENOME_TAXONOMY_QC } from '../subworkflows/local/genome_taxonomy_qc'
include { MOBILE_ELEMENTS }    from '../subworkflows/local/mobile_elements'

include { FASTP }                       from '../modules/nf-core/fastp/main'
include { SPADES }                      from '../modules/nf-core/spades/main'
include { CHECKM2_PREDICT }             from '../modules/nf-core/checkm2/predict/main'
include { PREP_ASSEMBLY }               from '../modules/local/util'
include { AVIARY_RECOVER }              from '../modules/local/aviary'
include { COVERM_CLUSTER; COVERM_GENOME } from '../modules/local/coverm'
include { PYRODIGAL as PYRODIGAL_SCAFFOLDS } from '../modules/local/pyrodigal'
include { NONPAREIL }                   from '../modules/local/nonpareil'

// resolve a possibly-null db param to a path or an empty list (for optional inputs)
def optpath = { p -> p ? file(p, checkIfExists: true) : [] }

workflow ILLUMINA_METAGENOME {

    if (!params.input) { error "Mode 'illumina_metagenome' requires --input <samplesheet.csv>" }

    INPUT_CHECK(params.input, params.mode)
    ch_reads = INPUT_CHECK.out.reads_short

    // --- QC ---
    if (!params.skip_qc) {
        FASTP(ch_reads.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
        ch_qc = FASTP.out.reads
    } else {
        ch_qc = ch_reads
    }

    // --- Read profiling (on raw reads, as per the bash workflow) ---
    ch_sylph_tax_meta = params.sylph_tax_metadata
        ? Channel.fromPath(params.sylph_tax_metadata, checkIfExists: true).collect()
        : Channel.value([])
    READ_PROFILING(
        ch_reads,
        Channel.fromPath(params.sylph_db, checkIfExists: true).collect(),
        ch_sylph_tax_meta,
        file(params.singlem_metapackage, checkIfExists: true),
        !params.skip_sylph,
        !params.skip_sylph && params.sylph_tax_metadata != null,
        !params.skip_singlem
    )

    // --- Host removal ---
    if (!params.skip_host_removal) {
        HOST_REMOVAL(ch_qc, optpath(params.host_ref), optpath(params.cleanifier_db))
        ch_clean = HOST_REMOVAL.out.reads
    } else {
        ch_clean = ch_qc
    }

    // --- Assembly (metaSPAdes) ---
    SPADES(ch_clean.map { meta, reads -> [ meta, reads, [], [] ] }, [], [])
    PREP_ASSEMBLY(SPADES.out.scaffolds)
    ch_assembly = PREP_ASSEMBLY.out.assembly        // [ meta, scaffolds.fasta ]

    // --- Binning (Aviary) ---
    ch_aviary_in = ch_assembly.join(ch_clean)       // [ meta, assembly, reads ]
    if (!params.skip_binning) {
        AVIARY_RECOVER(
            ch_aviary_in,
            file(params.gtdbtk_db,  checkIfExists: true),
            file(params.checkm2_db, checkIfExists: true),
            file(params.eggnog_db,  checkIfExists: true)
        )
        // all bins across samples, flattened for cross-sample steps
        ch_all_bins = AVIARY_RECOVER.out.bins.map { meta, bins -> bins }.flatten().collect()

        // CheckM2 on all bins (quality report drives dereplication)
        CHECKM2_PREDICT(
            ch_all_bins.map { bins -> [ [id: 'all_bins'], bins ] },
            [ [:], file(params.checkm2_db) ]
        )

        // --- Dereplication ---
        if (!params.skip_dereplication) {
            COVERM_CLUSTER(ch_all_bins, CHECKM2_PREDICT.out.checkm2_tsv.map { m, t -> t })
            ch_reps      = COVERM_CLUSTER.out.representatives.collect()
            ch_per_rep   = COVERM_CLUSTER.out.representatives.flatten()
                                .map { f -> [ [id: f.baseName], f ] }
        } else {
            ch_reps    = ch_all_bins
            ch_per_rep = AVIARY_RECOVER.out.bins.transpose().map { m, b -> [ [id: b.baseName], b ] }
        }

        // --- Map QC'd reads to representatives ---
        if (!params.skip_read_mapping) {
            COVERM_GENOME(ch_clean, ch_reps)
        }

        // --- Taxonomy + QC on representatives ---
        GENOME_TAXONOMY_QC(
            ch_reps,
            ch_per_rep,
            file(params.gtdbtk_db, checkIfExists: true),
            optpath(params.genomespot_models),
            !params.skip_taxonomy,
            params.run_checkm1,
            params.run_genomespot
        )
    }

    // --- Gene catalogue (from scaffold proteins) ---
    if (!params.skip_gene_catalogue) {
        PYRODIGAL_SCAFFOLDS(ch_assembly)
        GENE_CATALOGUE(
            PYRODIGAL_SCAFFOLDS.out.faa,
            PYRODIGAL_SCAFFOLDS.out.fna,
            optpath(params.dram_db),
            !params.skip_annotation
        )
    }

    // --- Mobile elements (virus/plasmid) ---
    if (!params.skip_mobile_elements) {
        MOBILE_ELEMENTS(
            ch_assembly,
            file(params.genomad_db, checkIfExists: true),
            file(params.checkv_db,  checkIfExists: true)
        )
    }

    // --- Sequencing coverage assessment ---
    if (params.run_nonpareil) {
        NONPAREIL(ch_clean)
    }
}
