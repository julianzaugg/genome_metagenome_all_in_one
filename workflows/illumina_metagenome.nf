/*
 * ILLUMINA_METAGENOME — full end-to-end workflow.
 *
 *  fastp QC ─┬─> sylph + singlem (read profiling, on raw reads)
 *            └─> host removal ─> metaSPAdes ─┬─> CoverM mapping (per-sample scaffolds)
 *                                            ├─> pyrodigal ─> gene catalogue ─> DRAM
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
include { AVIARY_RECOVER; AVIARY_COLLECT_BINS } from '../modules/local/aviary'
include { COVERM_CLUSTER; COVERM_GENOME; COVERM_CONTIG } from '../modules/local/coverm'
include { CHECKM1_LINEAGEWF }           from '../modules/local/checkm1'
include { PYRODIGAL as PYRODIGAL_SCAFFOLDS } from '../modules/local/pyrodigal'
include { NONPAREIL }                   from '../modules/local/nonpareil'
include { SEQKIT_STATS }                from '../modules/local/read_stats'
include { CLEANIFIER_INDEX }            from '../modules/local/host_removal'
include { FASTQ_GZIP_TEST }             from '../modules/local/validate'

// resolve a possibly-null db param to a path or an empty list (for optional inputs)
def optpath = { p -> p ? file(p, checkIfExists: true) : [] }

workflow ILLUMINA_METAGENOME {

    if (!params.input) { error "Mode 'illumina_metagenome' requires --input <samplesheet.csv>" }

    INPUT_CHECK(params.input, params.mode)
    FASTQ_GZIP_TEST(INPUT_CHECK.out.reads_short)
    ch_reads = FASTQ_GZIP_TEST.out.reads
    ch_read_stats = ch_reads.map { meta, reads -> [ meta, 'raw', reads ] }

    // --- QC ---
    if (!params.skip_qc) {
        FASTP(ch_reads.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
        ch_qc = FASTP.out.reads
        ch_read_stats = ch_read_stats.mix(FASTP.out.reads.map { meta, reads -> [ meta, 'fastp', reads ] })
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
        if (params.cleanifier_db) {
            ch_cleanifier_index = Channel.value(file(params.cleanifier_db, checkIfExists: true))
        } else {
            if (!params.host_ref) {
                error "Host removal is enabled but neither --cleanifier_db nor --host_ref is set. Provide a Cleanifier .filter index, provide a FASTA with --host_ref to build one, or rerun with --skip_host_removal true."
            }
            CLEANIFIER_INDEX(file(params.host_ref, checkIfExists: true), params.cleanifier_nobjects ?: '')
            ch_cleanifier_index = CLEANIFIER_INDEX.out.index
        }
        HOST_REMOVAL(ch_qc, ch_cleanifier_index)
        ch_clean = HOST_REMOVAL.out.reads
        ch_read_stats = ch_read_stats.mix(HOST_REMOVAL.out.reads.map { meta, reads -> [ meta, 'cleanifier', reads ] })
    } else {
        ch_clean = ch_qc
    }

    SEQKIT_STATS(ch_read_stats)

    // --- Assembly (metaSPAdes) ---
    ch_assembly = Channel.empty()
    if (!params.skip_assembly) {
        SPADES(ch_clean.map { meta, reads -> [ meta, reads, [], [] ] }, [], [])
        PREP_ASSEMBLY(SPADES.out.scaffolds)
        ch_assembly = PREP_ASSEMBLY.out.assembly    // [ meta, scaffolds.fasta ]
    }

    // --- Map QC'd reads to each sample's assembled scaffolds ---
    if (!params.skip_assembly && !params.skip_read_mapping) {
        COVERM_CONTIG(
            ch_assembly
                .join(ch_clean)
                .map { meta, scaffolds, reads -> [ meta, reads, scaffolds ] }
        )
    }

    // --- Binning (Aviary) ---
    ch_aviary_in = ch_assembly.join(ch_clean)       // [ meta, assembly, reads ]
    if (!params.skip_binning) {
        AVIARY_RECOVER(
            ch_aviary_in,
            file(params.gtdbtk_db,  checkIfExists: true),
            file(params.checkm2_db, checkIfExists: true),
            file(params.eggnog_db,  checkIfExists: true)
        )
        // all renamed bins across samples; this is the canonical bin set for
        // CheckM, dereplication, read mapping, and downstream taxonomy/QC.
        AVIARY_COLLECT_BINS(AVIARY_RECOVER.out.bins.map { meta, bins -> bins }.flatten().collect())
        ch_all_bins = AVIARY_COLLECT_BINS.out.bins

        // CheckM on ALL bins (pre-dereplication). CheckM2 drives clustering; both
        // CheckM1 and CheckM2 feed the high-quality-representative selection.
        ch_checkm2_tsv = Channel.value([])
        ch_checkm1_tsv = Channel.value([])
        if (!params.skip_checkm) {
            CHECKM2_PREDICT(
                ch_all_bins.map { bins -> [ [id: 'all_bins'], bins ] },
                [ [:], file(params.checkm2_db) ]
            )
            ch_checkm2_tsv = CHECKM2_PREDICT.out.checkm2_tsv.map { m, t -> t }
        }
        if (params.run_checkm1) {
            CHECKM1_LINEAGEWF(ch_all_bins)
            ch_checkm1_tsv = CHECKM1_LINEAGEWF.out.summary
        }

        // --- Dereplication ---
        if (!params.skip_dereplication) {
            COVERM_CLUSTER(ch_all_bins, ch_checkm2_tsv, ch_checkm1_tsv)
            ch_reps      = COVERM_CLUSTER.out.representatives.collect()
            ch_per_rep   = COVERM_CLUSTER.out.representatives.flatten()
                                .map { f -> [ [id: f.baseName], f ] }
        } else {
            ch_reps    = ch_all_bins
            ch_per_rep = AVIARY_COLLECT_BINS.out.bins.flatten().map { b -> [ [id: b.baseName], b ] }
        }

        // --- Map QC'd reads to representatives ---
        if (!params.skip_read_mapping) {
            COVERM_GENOME(ch_clean, ch_reps)
        }

        // --- Taxonomy + per-genome QC on representatives ---
        GENOME_TAXONOMY_QC(
            ch_reps,
            ch_per_rep,
            file(params.gtdbtk_db, checkIfExists: true),
            optpath(params.genomespot_models),
            !params.skip_taxonomy,
            params.run_genomespot,
            params.run_barrnap
        )
    }

    // --- Gene catalogue (from scaffold proteins) ---
    if (!params.skip_gene_catalogue) {
        PYRODIGAL_SCAFFOLDS(ch_assembly)
        GENE_CATALOGUE(
            PYRODIGAL_SCAFFOLDS.out.faa,
            PYRODIGAL_SCAFFOLDS.out.fna,
            params.catalogue_identities,
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
