/*
 * COMPARISON_READS — map an external set of reads (never assembled or
 * binned) against what this run recovered, for cross-cohort comparison:
 *   - sylph + singlem community profiling (raw reads, mirrors READ_PROFILING)
 *   - QC + host removal (mirrors the main sample path)
 *   - CoverM mapping against the final dereplicated bin representatives
 *   - RPKM-style DIAMOND mapping against the gene catalogue (base or
 *     expanded tier, chosen by the caller)
 *
 * Every process here is an ALIASED import of an already-existing process
 * (X as X_COMPARISON) — required, not stylistic: several of the underlying
 * processes emit collated, non-meta-keyed output (SYLPH_PROFILE,
 * SINGLEM_PIPE, RPKM_COLLATE_GENE_BLAST, RPKM_CALCULATE), so invoking the
 * un-aliased process from this second call site would publish over the
 * main run's own output of the same static filename.
 */

include { SEQKIT_STATS as SEQKIT_STATS_COMPARISON }         from '../../modules/local/read_stats'
include { SYLPH_SKETCH as SYLPH_SKETCH_COMPARISON;
          SYLPH_PROFILE as SYLPH_PROFILE_COMPARISON;
          SYLPH_TAX as SYLPH_TAX_COMPARISON }                from '../../modules/local/sylph'
include { SINGLEM_PIPE as SINGLEM_PIPE_COMPARISON }          from '../../modules/local/singlem'
include { FASTP as FASTP_COMPARISON }                        from '../../modules/nf-core/fastp/main'
include { CLEANIFIER as CLEANIFIER_COMPARISON }               from '../../modules/local/host_removal'
include { COVERM_GENOME as COVERM_GENOME_COMPARISON }        from '../../modules/local/coverm'
include { MAPPING_ASSESS as MAPPING_ASSESS_COMPARISON }      from '../../modules/local/mapping_assessment'
include {
    RPKM_FILTER_READS as RPKM_FILTER_READS_COMPARISON;
    RPKM_SINGLEM_MARKERS as RPKM_SINGLEM_MARKERS_COMPARISON;
    RPKM_PREBUILT_SINGLEM_MARKERS as RPKM_PREBUILT_SINGLEM_MARKERS_COMPARISON;
    RPKM_DIAMOND_MAKEDB as RPKM_DIAMOND_MAKEDB_COMPARISON;
    RPKM_DIAMOND_GENE as RPKM_DIAMOND_GENE_COMPARISON;
    RPKM_COLLATE_GENE_BLAST as RPKM_COLLATE_GENE_BLAST_COMPARISON;
    RPKM_DIAMOND_SINGLEM as RPKM_DIAMOND_SINGLEM_COMPARISON;
    RPKM_CALCULATE as RPKM_CALCULATE_COMPARISON
} from '../../modules/local/rpkm'

workflow COMPARISON_READS {
    take:
    reads                   // [ meta, [fastq_1, fastq_2] ]  from a second INPUT_CHECK invocation
    run_qc                  // bool = !params.skip_qc
    run_host_removal        // bool = !params.skip_host_removal
    cleanifier_index        // value channel [filter, info]; ignored if run_host_removal is false
    run_sylph
    run_sylph_tax
    sylph_db
    sylph_tax_metadata
    run_singlem
    singlem_metapackage
    run_bin_mapping          // bool = !params.skip_dereplication && !params.skip_read_mapping
    bin_representatives      // ch_reps (collected fasta paths)
    run_mapping_assessment   // bool = !params.skip_mapping_assessment
    run_gene_mapping         // bool = !params.skip_gene_catalogue
    gene_catalogue           // base or expanded GENE_CATALOGUE output, chosen by the caller
    singlem_marker_dbs
    singlem_marker_lengths
    min_read_length

    main:
    ch_versions = Channel.empty()

    // Raw-read totals -- the denominator MAPPING_ASSESS needs.
    SEQKIT_STATS_COMPARISON(reads.map { meta, r -> [ meta, 'raw', r ] })
    ch_raw_stats = SEQKIT_STATS_COMPARISON.out.stats.map { meta, stage, t -> [ meta, t ] }
    ch_versions = ch_versions.mix(SEQKIT_STATS_COMPARISON.out.versions)

    // --- Community profiling (raw reads; mirrors READ_PROFILING) ---
    ch_sylph_profile = Channel.empty()
    if (run_sylph) {
        SYLPH_SKETCH_COMPARISON(reads)
        ch_sketches = SYLPH_SKETCH_COMPARISON.out.sketch.map { meta, s -> s }.collect()
        SYLPH_PROFILE_COMPARISON(ch_sketches, sylph_db)
        ch_sylph_profile = SYLPH_PROFILE_COMPARISON.out.profile
        ch_versions = ch_versions
            .mix(SYLPH_SKETCH_COMPARISON.out.versions)
            .mix(SYLPH_PROFILE_COMPARISON.out.versions)

        if (run_sylph_tax) {
            SYLPH_TAX_COMPARISON(SYLPH_PROFILE_COMPARISON.out.profile, sylph_tax_metadata)
            ch_versions = ch_versions.mix(SYLPH_TAX_COMPARISON.out.versions)
        }
    }

    ch_singlem_profile = Channel.empty()
    if (run_singlem) {
        ch_read_sets = reads.collect(flat: false)
        ch_singlem_forward = ch_read_sets.map { rows ->
            rows.collect { row ->
                def meta = row[0]
                def r    = row[1]
                meta.single_end ? r : r[0]
            }
        }
        ch_singlem_reverse = ch_read_sets.map { rows ->
            rows.findAll { row -> !row[0].single_end }.collect { row -> row[1][1] }
        }
        SINGLEM_PIPE_COMPARISON(ch_singlem_forward, ch_singlem_reverse, singlem_metapackage)
        ch_singlem_profile = SINGLEM_PIPE_COMPARISON.out.profile
        ch_versions = ch_versions.mix(SINGLEM_PIPE_COMPARISON.out.versions)
    }

    // --- QC ---
    if (run_qc) {
        FASTP_COMPARISON(reads.map { meta, r -> [ meta, r, [] ] }, false, false, false)
        ch_qc = FASTP_COMPARISON.out.reads
        // FASTP emits versions via topic: versions (collected globally in illumina_metagenome.nf)
    } else {
        ch_qc = reads
    }

    // --- Host removal ---
    if (run_host_removal) {
        CLEANIFIER_COMPARISON(ch_qc, cleanifier_index)
        ch_clean = CLEANIFIER_COMPARISON.out.reads
        ch_versions = ch_versions.mix(CLEANIFIER_COMPARISON.out.versions)
    } else {
        ch_clean = ch_qc
    }

    // --- Map to the final dereplicated bin representatives ---
    ch_bin_abundance = Channel.empty()
    if (run_bin_mapping) {
        COVERM_GENOME_COMPARISON(ch_clean, bin_representatives)
        ch_bin_abundance = COVERM_GENOME_COMPARISON.out.abundance
        ch_versions = ch_versions.mix(COVERM_GENOME_COMPARISON.out.versions)

        if (run_mapping_assessment) {
            MAPPING_ASSESS_COMPARISON(
                COVERM_GENOME_COMPARISON.out.bams.join(COVERM_GENOME_COMPARISON.out.contig_map).join(ch_raw_stats),
                'Comparison_vs_Bins'
            )
            ch_versions = ch_versions.mix(MAPPING_ASSESS_COMPARISON.out.versions)
        }
    }

    // --- Map to the gene catalogue (RPKM-style; own diamond db + own SingleM
    // marker blast, since these are different reads from the main RPKM run) ---
    ch_gene_rpkm = Channel.empty()
    if (run_gene_mapping) {
        ch_rpkm_r1 = ch_clean.map { meta, r -> [ meta, meta.single_end ? r : r[0] ] }
        RPKM_FILTER_READS_COMPARISON(ch_rpkm_r1, min_read_length)
        ch_versions = ch_versions.mix(RPKM_FILTER_READS_COMPARISON.out.versions)

        if (singlem_marker_dbs && singlem_marker_lengths) {
            RPKM_PREBUILT_SINGLEM_MARKERS_COMPARISON(singlem_marker_dbs, singlem_marker_lengths)
            ch_marker_dbs     = RPKM_PREBUILT_SINGLEM_MARKERS_COMPARISON.out.marker_dbs
            ch_marker_lengths = RPKM_PREBUILT_SINGLEM_MARKERS_COMPARISON.out.marker_lengths
            ch_versions = ch_versions.mix(RPKM_PREBUILT_SINGLEM_MARKERS_COMPARISON.out.versions)
        } else if (!singlem_marker_dbs && !singlem_marker_lengths) {
            RPKM_SINGLEM_MARKERS_COMPARISON(singlem_metapackage)
            ch_marker_dbs     = RPKM_SINGLEM_MARKERS_COMPARISON.out.marker_dbs
            ch_marker_lengths = RPKM_SINGLEM_MARKERS_COMPARISON.out.marker_lengths
            ch_versions = ch_versions.mix(RPKM_SINGLEM_MARKERS_COMPARISON.out.versions)
        } else {
            error "COMPARISON_READS requires both --rpkm_singlem_marker_dbs and --rpkm_singlem_marker_lengths when using prebuilt SingleM marker inputs."
        }

        RPKM_DIAMOND_MAKEDB_COMPARISON(gene_catalogue)
        ch_versions = ch_versions.mix(RPKM_DIAMOND_MAKEDB_COMPARISON.out.versions)

        ch_reads_for_gene = RPKM_FILTER_READS_COMPARISON.out.reads
            .combine(RPKM_DIAMOND_MAKEDB_COMPARISON.out.db)
            .map { meta, read, db -> [ meta, read, db ] }
        RPKM_DIAMOND_GENE_COMPARISON(ch_reads_for_gene)
        RPKM_COLLATE_GENE_BLAST_COMPARISON(RPKM_DIAMOND_GENE_COMPARISON.out.blast.map { meta, blast -> blast }.collect())
        ch_versions = ch_versions
            .mix(RPKM_DIAMOND_GENE_COMPARISON.out.versions)
            .mix(RPKM_COLLATE_GENE_BLAST_COMPARISON.out.versions)

        ch_reads_for_singlem = RPKM_FILTER_READS_COMPARISON.out.reads
            .combine(ch_marker_dbs)
            .map { meta, read, marker_dbs -> [ meta, read, marker_dbs ] }
        RPKM_DIAMOND_SINGLEM_COMPARISON(ch_reads_for_singlem)
        ch_versions = ch_versions.mix(RPKM_DIAMOND_SINGLEM_COMPARISON.out.versions)

        RPKM_CALCULATE_COMPARISON(
            RPKM_DIAMOND_SINGLEM_COMPARISON.out.blast_dir.collect(),
            ch_marker_lengths,
            RPKM_COLLATE_GENE_BLAST_COMPARISON.out.blast
        )
        ch_gene_rpkm = RPKM_CALCULATE_COMPARISON.out.normalised_rpkm
        ch_versions = ch_versions.mix(RPKM_CALCULATE_COMPARISON.out.versions)
    }

    emit:
    sylph_profile   = ch_sylph_profile
    singlem_profile = ch_singlem_profile
    bin_abundance   = ch_bin_abundance
    gene_rpkm       = ch_gene_rpkm
    versions        = ch_versions
}
