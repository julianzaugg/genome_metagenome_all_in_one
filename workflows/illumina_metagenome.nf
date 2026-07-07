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
include { REFERENCE_GENOMES }  from '../subworkflows/local/reference_genomes'
include { GENE_CATALOGUE }     from '../subworkflows/local/gene_catalogue'
include { GENOME_TAXONOMY_QC } from '../subworkflows/local/genome_taxonomy_qc'
include { MARKER_GENE_TREE }   from '../subworkflows/local/marker_gene_tree'
include { MOBILE_ELEMENTS }    from '../subworkflows/local/mobile_elements'
include { RPKM }               from '../subworkflows/local/rpkm'

include { FASTP }                       from '../modules/nf-core/fastp/main'
include { SPADES }                      from '../modules/nf-core/spades/main'
include { CHECKM2_PREDICT }             from '../modules/nf-core/checkm2/predict/main'
include { PREP_ASSEMBLY }               from '../modules/local/util'
include { AVIARY_RECOVER; AVIARY_COLLECT_BINS } from '../modules/local/aviary'
include { COVERM_CLUSTER; COVERM_CLUSTER_HQ; COVERM_CLUSTER_HQ_REF; COVERM_GENOME; COVERM_GENOME as COVERM_GENOME_HQ; COVERM_GENOME as COVERM_GENOME_HQ_DEREP; COVERM_GENOME as COVERM_GENOME_HQ_REF; COVERM_CONTIG } from '../modules/local/coverm'
include { CHECKM1_LINEAGEWF }           from '../modules/local/checkm1'
include { PYRODIGAL as PYRODIGAL_SCAFFOLDS } from '../modules/local/pyrodigal'
include { NONPAREIL }                   from '../modules/local/nonpareil'
include { SEQKIT_STATS }                from '../modules/local/read_stats'
include { READ_STAT_REPORT }            from '../modules/local/read_stat_report'
include { CLEANIFIER_INDEX }            from '../modules/local/host_removal'
include { FASTQ_GZIP_TEST }             from '../modules/local/validate'
include { DUMP_SOFTWARE_VERSIONS }      from '../modules/local/dump_software_versions'

// resolve a possibly-null db param to a path or an empty list (for optional inputs)
def optpath = { p -> p ? file(p, checkIfExists: true) : [] }

workflow ILLUMINA_METAGENOME {

    if (!params.input) { error "Mode 'illumina_metagenome' requires --input <samplesheet.csv>" }
    if (!params.skip_rpkm && params.skip_qc) {
        error "RPKM is enabled but --skip_qc true would bypass fastp. Run fastp or rerun with --skip_rpkm true."
    }
    if (!params.skip_rpkm && (params.skip_assembly || params.skip_gene_catalogue)) {
        error "RPKM requires assembly and the gene catalogue. Rerun with --skip_assembly false --skip_gene_catalogue false, or set --skip_rpkm true."
    }
    if (params.reference_genomes) {
        if (params.skip_binning || params.skip_dereplication) {
            error "--reference_genomes dereplicates the references with the HQ MAGs, which needs binning + dereplication. Rerun with --skip_binning false --skip_dereplication false, or unset --reference_genomes."
        }
        if (params.skip_checkm) {
            error "--reference_genomes selects cluster representatives by CheckM2 quality, so CheckM2 is required. Rerun with --skip_checkm false, or unset --reference_genomes."
        }
    }

    INPUT_CHECK(params.input, params.mode)
    FASTQ_GZIP_TEST(INPUT_CHECK.out.reads_short)
    ch_reads = FASTQ_GZIP_TEST.out.reads
    ch_read_stats = ch_reads.map { meta, reads -> [ meta, 'raw', reads ] }
    ch_versions = Channel.empty()
        .mix(INPUT_CHECK.out.versions)
        .mix(FASTQ_GZIP_TEST.out.versions)

    // --- QC ---
    if (!params.skip_qc) {
        FASTP(ch_reads.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
        ch_qc = FASTP.out.reads
        ch_read_stats = ch_read_stats.mix(FASTP.out.reads.map { meta, reads -> [ meta, 'fastp', reads ] })
        // FASTP emits versions via topic: versions (collected globally below)
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
    ch_versions = ch_versions.mix(READ_PROFILING.out.versions)

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
            ch_versions = ch_versions.mix(CLEANIFIER_INDEX.out.versions)
        }
        HOST_REMOVAL(ch_qc, ch_cleanifier_index)
        ch_clean = HOST_REMOVAL.out.reads
        ch_read_stats = ch_read_stats.mix(HOST_REMOVAL.out.reads.map { meta, reads -> [ meta, 'cleanifier', reads ] })
        ch_versions = ch_versions.mix(HOST_REMOVAL.out.versions)
    } else {
        ch_clean = ch_qc
    }

    SEQKIT_STATS(ch_read_stats)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions)

    // Read-stat report inputs (filled in as the relevant steps run)
    ch_scaffold_counts = Channel.value([])
    ch_repmag_abund    = Channel.value([])
    ch_hq_reps         = Channel.value([])
    ch_hq_repmag_abund = Channel.value([])
    ch_hq_derep_abund  = Channel.value([])
    ch_hq_ref_abund    = Channel.value([])

    // --- External reference genomes (normalise + CheckM2 + protein prediction) ---
    if (params.reference_genomes) {
        ch_reference_files = Channel
            .fromPath("${params.reference_genomes}/*.{${params.reference_genome_extension}}", checkIfExists: true)
            .collect()
        REFERENCE_GENOMES(
            ch_reference_files,
            optpath(params.reference_genomes_checkm2),
            file(params.checkm2_db, checkIfExists: true),
            params.reference_genomes_checkm2 == null
        )
        ch_versions = ch_versions.mix(REFERENCE_GENOMES.out.versions)
    }

    // --- Assembly (metaSPAdes) ---
    ch_assembly = Channel.empty()
    if (!params.skip_assembly) {
        SPADES(ch_clean.map { meta, reads -> [ meta, reads, [], [] ] }, [], [])
        PREP_ASSEMBLY(SPADES.out.scaffolds)
        ch_assembly = PREP_ASSEMBLY.out.assembly    // [ meta, scaffolds.fasta ]
        // SPADES emits versions via topic: versions (collected globally below)
    }

    // --- Map QC'd reads to each sample's assembled scaffolds ---
    if (!params.skip_assembly && !params.skip_read_mapping) {
        COVERM_CONTIG(
            ch_assembly
                .join(ch_clean)
                .map { meta, scaffolds, reads -> [ meta, reads, scaffolds ] }
        )
        ch_scaffold_counts = COVERM_CONTIG.out.counts.map { meta, t -> t }.collect().ifEmpty([])
        ch_versions = ch_versions.mix(COVERM_CONTIG.out.versions)
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
        ch_versions = ch_versions
            .mix(AVIARY_RECOVER.out.versions)
            .mix(AVIARY_COLLECT_BINS.out.versions)

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
            // CHECKM2_PREDICT emits versions via topic: versions (collected globally below)
        }
        if (params.run_checkm1) {
            CHECKM1_LINEAGEWF(ch_all_bins)
            ch_checkm1_tsv = CHECKM1_LINEAGEWF.out.summary
            ch_versions = ch_versions.mix(CHECKM1_LINEAGEWF.out.versions)
        }

        // --- Dereplication ---
        if (!params.skip_dereplication) {
            COVERM_CLUSTER(ch_all_bins, ch_checkm2_tsv, ch_checkm1_tsv)
            ch_reps      = COVERM_CLUSTER.out.representatives.collect()
            ch_per_rep   = COVERM_CLUSTER.out.representatives.flatten()
                                .map { f -> [ [id: f.baseName], f ] }
            ch_hq_reps   = COVERM_CLUSTER.out.hq_representatives.collect().ifEmpty([])
            ch_versions = ch_versions.mix(COVERM_CLUSTER.out.versions)

            // HQ-first: extract HQ MAGs from the FULL bin set, THEN dereplicate them
            COVERM_CLUSTER_HQ(ch_all_bins, ch_checkm2_tsv, ch_checkm1_tsv)
            ch_hq_derep_reps = COVERM_CLUSTER_HQ.out.representatives.collect().ifEmpty([])
            ch_versions = ch_versions.mix(COVERM_CLUSTER_HQ.out.versions)

            // HQ-first + external reference genomes: dereplicate the HQ MAGs together with
            // the references (references always included; representatives chosen by CheckM2).
            if (params.reference_genomes) {
                COVERM_CLUSTER_HQ_REF(ch_all_bins, REFERENCE_GENOMES.out.genomes,
                                      ch_checkm2_tsv, ch_checkm1_tsv, REFERENCE_GENOMES.out.checkm2)
                ch_hq_ref_derep_reps = COVERM_CLUSTER_HQ_REF.out.representatives.collect().ifEmpty([])
                ch_versions = ch_versions.mix(COVERM_CLUSTER_HQ_REF.out.versions)
            }
        } else {
            ch_reps    = ch_all_bins
            ch_per_rep = AVIARY_COLLECT_BINS.out.bins.flatten().map { b -> [ [id: b.baseName], b ] }
        }

        // --- Map QC'd reads to representatives ---
        if (!params.skip_read_mapping) {
            COVERM_GENOME(ch_clean, ch_reps)
            ch_repmag_abund = COVERM_GENOME.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
            ch_versions = ch_versions.mix(COVERM_GENOME.out.versions)

            // --- Map the same reads directly to the HQ-only subset (not extracted from
            // the full-set mapping above) — avoids undercounting HQ MAGs when dereplication
            // leaves redundant lower-quality near-duplicate bins that split reads away from them ---
            if (!params.skip_dereplication) {
                COVERM_GENOME_HQ(ch_clean, ch_hq_reps)
                ch_hq_repmag_abund = COVERM_GENOME_HQ.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
                ch_versions = ch_versions.mix(COVERM_GENOME_HQ.out.versions)

                // Map the same reads to the HQ-first-then-dereplicated set
                COVERM_GENOME_HQ_DEREP(ch_clean, ch_hq_derep_reps)
                ch_hq_derep_abund = COVERM_GENOME_HQ_DEREP.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
                ch_versions = ch_versions.mix(COVERM_GENOME_HQ_DEREP.out.versions)

                // Map the same reads to the (HQ MAGs + reference genomes) dereplicated set
                if (params.reference_genomes) {
                    COVERM_GENOME_HQ_REF(ch_clean, ch_hq_ref_derep_reps)
                    ch_hq_ref_abund = COVERM_GENOME_HQ_REF.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
                    ch_versions = ch_versions.mix(COVERM_GENOME_HQ_REF.out.versions)
                }
            }
        }

        // --- Taxonomy + per-genome QC on all bins ---
        GENOME_TAXONOMY_QC(
            ch_all_bins,
            AVIARY_COLLECT_BINS.out.bins.flatten().map { b -> [ [id: b.baseName], b ] },
            file(params.gtdbtk_db, checkIfExists: true),
            optpath(params.genomespot_models),
            optpath(params.dram_db),
            !params.skip_taxonomy,
            params.run_genomespot,
            params.run_barrnap,
            params.run_dram_bins
        )
        ch_versions = ch_versions.mix(GENOME_TAXONOMY_QC.out.versions)

        // --- Marker-gene tree of MAGs + selected GTDB references ---
        if (params.run_marker_tree && !params.skip_taxonomy) {
            MARKER_GENE_TREE(
                GENOME_TAXONOMY_QC.out.gtdbtk_outdir,
                params.marker_tree_genome_source == 'all_bins' ? ch_all_bins : ch_reps,
                ch_checkm2_tsv
            )
            ch_versions = ch_versions.mix(MARKER_GENE_TREE.out.versions)
        }
    }

    // --- Gene catalogue (from scaffold proteins) ---
    if (!params.skip_gene_catalogue) {
        PYRODIGAL_SCAFFOLDS(ch_assembly)
        GENE_CATALOGUE(
            PYRODIGAL_SCAFFOLDS.out.faa,
            PYRODIGAL_SCAFFOLDS.out.fna,
            params.catalogue_identities,
            optpath(params.dram_db),
            !params.skip_annotation,
            params.reference_genomes ? REFERENCE_GENOMES.out.faa : Channel.empty(),
            params.reference_genomes ? REFERENCE_GENOMES.out.fna : Channel.empty(),
            params.reference_genomes as boolean
        )
        ch_versions = ch_versions
            .mix(PYRODIGAL_SCAFFOLDS.out.versions)
            .mix(GENE_CATALOGUE.out.versions)
    }

    // --- RPKM (selected stream only: host-filtered fastp reads, or fastp reads if host removal is skipped) ---
    if (!params.skip_rpkm) {
        ch_rpkm_r1 = ch_clean.map { meta, reads -> [ meta, meta.single_end ? reads : reads[0] ] }
        RPKM(
            ch_rpkm_r1,
            GENE_CATALOGUE.out.catalogue,
            file(params.singlem_metapackage, checkIfExists: true),
            optpath(params.rpkm_singlem_marker_dbs),
            optpath(params.rpkm_singlem_marker_lengths),
            params.rpkm_min_read_length,
            params.reference_genomes ? GENE_CATALOGUE.out.expanded_catalogue : []
        )
        ch_versions = ch_versions.mix(RPKM.out.versions)
    }

    // --- Mobile elements (virus/plasmid) ---
    if (!params.skip_mobile_elements) {
        MOBILE_ELEMENTS(
            ch_assembly,
            file(params.genomad_db, checkIfExists: true),
            file(params.checkv_db,  checkIfExists: true)
        )
        ch_versions = ch_versions.mix(MOBILE_ELEMENTS.out.versions)
    }

    // --- Sequencing coverage assessment ---
    if (params.run_nonpareil) {
        NONPAREIL(ch_clean)
        ch_versions = ch_versions.mix(NONPAREIL.out.versions)
    }

    // --- Read-stat report (per-sample read tracking across all steps) ---
    READ_STAT_REPORT(
        'metagenome',
        SEQKIT_STATS.out.stats.map { meta, stage, t -> t }.collect(),
        ch_scaffold_counts,
        [],
        ch_repmag_abund,
        ch_hq_reps,
        ch_hq_repmag_abund,
        ch_hq_derep_abund,
        ch_hq_ref_abund
    )
    ch_versions = ch_versions.mix(READ_STAT_REPORT.out.versions)

    // --- Software versions manifest ---
    // nf-core modules emit tuples via topic: versions; convert to file and merge with local versions
    ch_nfcore_versions = channel.topic('versions')
        .collectFile(name: 'nfcore_versions.yml', newLine: true) { process, tool, version ->
            "\"${process}\":\n    ${tool}: ${version}"
        }
    DUMP_SOFTWARE_VERSIONS(ch_versions.mix(ch_nfcore_versions).collect())
}
