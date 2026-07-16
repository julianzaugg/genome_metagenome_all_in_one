/*
 * NANOPORE_METAGENOME — long-read metagenome workflow.
 */

include { INPUT_CHECK }        from '../subworkflows/local/input_check'
include { LONG_READ_QC }       from '../subworkflows/local/long_read_qc'
include { READ_PROFILING }     from '../subworkflows/local/read_profiling'
include { HOST_REMOVAL }       from '../subworkflows/local/host_removal'
include { REFERENCE_GENOMES }  from '../subworkflows/local/reference_genomes'
include { GENE_CATALOGUE }     from '../subworkflows/local/gene_catalogue'
include { GENOME_TAXONOMY_QC } from '../subworkflows/local/genome_taxonomy_qc'
include { MARKER_GENE_TREE }   from '../subworkflows/local/marker_gene_tree'
include { MOBILE_ELEMENTS }    from '../subworkflows/local/mobile_elements'

include { MYLOASM }                     from '../modules/local/assembly_isolate'
include { DORADO_POLISH }               from '../modules/local/long_reads'
include { CHECKM2_PREDICT }             from '../modules/nf-core/checkm2/predict/main'
include { AVIARY_RECOVER; AVIARY_COLLECT_BINS } from '../modules/local/aviary'
include { COVERM_CLUSTER; COVERM_CLUSTER_HQ; COVERM_CLUSTER_HQ_REF; COVERM_GENOME as COVERM_GENOME_ONT; COVERM_GENOME as COVERM_GENOME_HQ_ONT; COVERM_GENOME as COVERM_GENOME_HQ_DEREP_ONT; COVERM_GENOME as COVERM_GENOME_HQ_REF_ONT; COVERM_CONTIG as COVERM_CONTIG_ONT } from '../modules/local/coverm'
include { COVERM_CLUSTER_WS; COVERM_CLUSTER_HQ_WS; COVERM_GENOME_PAIRED as COVERM_GENOME_WS_DEREP_ONT; COVERM_GENOME_PAIRED as COVERM_GENOME_WS_HQ_ONT } from '../modules/local/coverm'
include { CHECKM1_LINEAGEWF }           from '../modules/local/checkm1'
include { PYRODIGAL as PYRODIGAL_SCAFFOLDS } from '../modules/local/pyrodigal'
include { NONPAREIL }                   from '../modules/local/nonpareil'
include { READ_STAT_REPORT }            from '../modules/local/read_stat_report'
include { CLEANIFIER_INDEX }            from '../modules/local/host_removal'
include { DUMP_SOFTWARE_VERSIONS }      from '../modules/local/dump_software_versions'

def optpath = { p -> p ? file(p, checkIfExists: true) : [] }

workflow NANOPORE_METAGENOME {

    if (!params.input) { error "Mode 'nanopore_metagenome' requires --input <samplesheet.csv>" }
    if (params.reference_genomes) {
        if (params.skip_binning || params.skip_dereplication) {
            error "--reference_genomes dereplicates the references with the HQ MAGs, which needs binning + dereplication. Rerun with --skip_binning false --skip_dereplication false, or unset --reference_genomes."
        }
        if (params.skip_checkm) {
            error "--reference_genomes selects cluster representatives by CheckM2 quality, so CheckM2 is required. Rerun with --skip_checkm false, or unset --reference_genomes."
        }
    }
    if (params.reference_genomes_taxonomy && (!params.reference_genomes || params.skip_taxonomy)) {
        error "--reference_genomes_taxonomy needs --reference_genomes and taxonomy enabled (not --skip_taxonomy)."
    }
    if (params.marker_tree_include_references && (!params.reference_genomes || !params.run_marker_tree || params.skip_taxonomy)) {
        error "--marker_tree_include_references needs --reference_genomes, --run_marker_tree true, and taxonomy enabled (not --skip_taxonomy)."
    }
    if (params.marker_tree_genome_source == 'hq_representatives' && (params.skip_binning || params.skip_dereplication)) {
        error "--marker_tree_genome_source hq_representatives uses the high-quality representative set, which needs binning + dereplication. Rerun with --skip_binning false --skip_dereplication false, or choose another source."
    }
    if (params.within_sample_dereplication != 'none') {
        if (params.skip_binning) {
            error "--within_sample_dereplication needs binning. Rerun with --skip_binning false, or set --within_sample_dereplication none."
        }
        if (params.skip_checkm && !params.run_checkm1) {
            log.warn "within_sample_dereplication: no CheckM report (skip_checkm true, run_checkm1 false) -> within-sample HQ MAG selection will be empty."
        }
    }

    INPUT_CHECK(params.input, params.mode)

    ch_direct_long = INPUT_CHECK.out.reads_long
        .filter { meta, reads -> !params.force_dorado_basecalling || !meta.has_pod5 }
    ch_pod5_basecall = INPUT_CHECK.out.pod5
        .filter { meta, pod5 -> params.force_dorado_basecalling || !meta.has_long_reads }

    LONG_READ_QC(
        ch_direct_long,
        ch_pod5_basecall,
        params.dorado_model,
        params.dorado_barcode_kit ?: '',
        params.dorado_device,
        true,
        !params.skip_porechop,
        !params.skip_qc
    )
    ch_reads = LONG_READ_QC.out.reads
    ch_versions = Channel.empty()
        .mix(INPUT_CHECK.out.versions)
        .mix(LONG_READ_QC.out.versions)

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
        HOST_REMOVAL(ch_reads, ch_cleanifier_index)
        ch_clean = HOST_REMOVAL.out.reads
        ch_versions = ch_versions.mix(HOST_REMOVAL.out.versions)
    } else {
        ch_clean = ch_reads
    }

    // Read-stat report inputs (filled in as the relevant steps run)
    ch_scaffold_counts = Channel.value([])
    ch_repmag_abund    = Channel.value([])
    ch_hq_reps         = Channel.value([])
    ch_hq_repmag_abund = Channel.value([])
    ch_hq_derep_abund  = Channel.value([])
    ch_hq_ref_abund    = Channel.value([])
    ch_ws_derep_abund  = Channel.value([])
    ch_ws_hq_abund     = Channel.value([])

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

    ch_assembly = Channel.empty()
    if (!params.skip_assembly) {
        MYLOASM(ch_clean)
        ch_assembly = MYLOASM.out.assembly
        ch_versions = ch_versions.mix(MYLOASM.out.versions)

        if (!params.skip_dorado_polish) {
            DORADO_POLISH(
                ch_assembly.join(ch_clean).map { meta, assembly, reads -> [ meta, assembly, reads ] },
                params.dorado_device
            )
            ch_assembly = DORADO_POLISH.out.assembly
            ch_versions = ch_versions.mix(DORADO_POLISH.out.versions)
        }
    }

    if (!params.skip_assembly && !params.skip_read_mapping) {
        COVERM_CONTIG_ONT(
            ch_assembly
                .join(ch_clean)
                .map { meta, scaffolds, reads -> [ meta, reads, scaffolds ] }
        )
        ch_scaffold_counts = COVERM_CONTIG_ONT.out.counts.map { meta, t -> t }.collect().ifEmpty([])
        ch_versions = ch_versions.mix(COVERM_CONTIG_ONT.out.versions)
    }

    if (!params.skip_binning) {
        AVIARY_RECOVER(
            ch_assembly.join(ch_clean),
            file(params.gtdbtk_db,  checkIfExists: true),
            file(params.checkm2_db, checkIfExists: true),
            file(params.eggnog_db,  checkIfExists: true)
        )
        AVIARY_COLLECT_BINS(AVIARY_RECOVER.out.bins.map { meta, bins -> bins }.flatten().collect())
        ch_all_bins = AVIARY_COLLECT_BINS.out.bins
        ch_versions = ch_versions
            .mix(AVIARY_RECOVER.out.versions)
            .mix(AVIARY_COLLECT_BINS.out.versions)

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

        // --- Within-sample / within-group dereplication (independent of the cross-sample
        // path below; reuses the pooled CheckM reports — only each unit's bins are clustered) ---
        if (params.within_sample_dereplication != 'none') {
            def ck2 = ch_checkm2_tsv.first()   // broadcast the single-emission CheckM reports
            def ck1 = ch_checkm1_tsv.first()

            // grouping unit -> [meta, bins]; 'sample' keeps FULL meta (needed for the read join)
            ch_ws_bins = params.within_sample_dereplication == 'group'
                ? AVIARY_RECOVER.out.bins.map { meta, bins -> [ meta.group, bins ] }
                      .groupTuple().map { g, bl -> [ [id: g], bl.flatten() ] }
                : AVIARY_RECOVER.out.bins

            COVERM_CLUSTER_WS(ch_ws_bins, ck2, ck1)
            COVERM_CLUSTER_HQ_WS(ch_ws_bins, ck2, ck1)
            ch_versions = ch_versions.mix(COVERM_CLUSTER_WS.out.versions, COVERM_CLUSTER_HQ_WS.out.versions)

            // each sample maps ITS reads to its own unit's reps
            def joinReps = { repsCh ->
                params.within_sample_dereplication == 'group'
                    ? ch_clean.map { meta, reads -> [ meta.group, meta, reads ] }
                          .combine(repsCh.map { m, r -> [ m.id, r ] }, by: 0)
                          .map { grp, meta, reads, r -> [ meta, reads, r ] }
                    : ch_clean.join(repsCh)   // [meta, reads, reps]
            }

            if (!params.skip_read_mapping) {
                COVERM_GENOME_WS_DEREP_ONT(joinReps(COVERM_CLUSTER_WS.out.representatives))
                COVERM_GENOME_WS_HQ_ONT(joinReps(COVERM_CLUSTER_HQ_WS.out.representatives))
                ch_ws_derep_abund = COVERM_GENOME_WS_DEREP_ONT.out.abundance.map { m, t -> t }.collect().ifEmpty([])
                ch_ws_hq_abund    = COVERM_GENOME_WS_HQ_ONT.out.abundance.map { m, t -> t }.collect().ifEmpty([])
                ch_versions = ch_versions.mix(COVERM_GENOME_WS_DEREP_ONT.out.versions, COVERM_GENOME_WS_HQ_ONT.out.versions)
            }
        }

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

            // HQ-first + external reference genomes (references always included;
            // representatives chosen by CheckM2). See illumina_metagenome.nf.
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

        if (!params.skip_read_mapping) {
            COVERM_GENOME_ONT(ch_clean, ch_reps)
            ch_repmag_abund = COVERM_GENOME_ONT.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
            ch_versions = ch_versions.mix(COVERM_GENOME_ONT.out.versions)

            // --- Map the same reads directly to the HQ-only subset (see illumina_metagenome.nf
            // for why this is kept separate from the full-set subset extraction) ---
            if (!params.skip_dereplication) {
                COVERM_GENOME_HQ_ONT(ch_clean, ch_hq_reps)
                ch_hq_repmag_abund = COVERM_GENOME_HQ_ONT.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
                ch_versions = ch_versions.mix(COVERM_GENOME_HQ_ONT.out.versions)

                // Map the same reads to the HQ-first-then-dereplicated set
                COVERM_GENOME_HQ_DEREP_ONT(ch_clean, ch_hq_derep_reps)
                ch_hq_derep_abund = COVERM_GENOME_HQ_DEREP_ONT.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
                ch_versions = ch_versions.mix(COVERM_GENOME_HQ_DEREP_ONT.out.versions)

                // Map the same reads to the (HQ MAGs + reference genomes) dereplicated set
                if (params.reference_genomes) {
                    COVERM_GENOME_HQ_REF_ONT(ch_clean, ch_hq_ref_derep_reps)
                    ch_hq_ref_abund = COVERM_GENOME_HQ_REF_ONT.out.abundance.map { meta, t -> t }.collect().ifEmpty([])
                    ch_versions = ch_versions.mix(COVERM_GENOME_HQ_REF_ONT.out.versions)
                }
            }
        }

        // Classify the external reference genomes in the same GTDB-Tk run when the user
        // wants their taxonomy and/or wants them placed in the marker tree.
        def classify_refs = params.reference_genomes && (params.reference_genomes_taxonomy || params.marker_tree_include_references)
        GENOME_TAXONOMY_QC(
            ch_all_bins,
            AVIARY_COLLECT_BINS.out.bins.flatten().map { b -> [ [id: b.baseName], b ] },
            file(params.gtdbtk_db, checkIfExists: true),
            optpath(params.genomespot_models),
            optpath(params.dram_db),
            !params.skip_taxonomy,
            params.run_genomespot,
            params.run_barrnap,
            params.run_dram_bins,
            classify_refs ? REFERENCE_GENOMES.out.gtdbtk_genomes : Channel.value([]),
            classify_refs
        )
        ch_versions = ch_versions.mix(GENOME_TAXONOMY_QC.out.versions)

        // --- Marker-gene tree of MAGs + selected GTDB references ---
        if (params.run_marker_tree && !params.skip_taxonomy) {
            def mt_genomes = params.marker_tree_genome_source == 'all_bins'          ? ch_all_bins
                           : params.marker_tree_genome_source == 'hq_representatives' ? ch_hq_reps
                           :                                                            ch_reps
            // hq_representatives is already HQ-filtered upstream (comp - 3*cont >= 50);
            // pass no CheckM2 report so the 90/5 threshold filter is not re-applied.
            def mt_checkm2 = params.marker_tree_genome_source == 'hq_representatives' ? Channel.value([]) : ch_checkm2_tsv
            MARKER_GENE_TREE(
                GENOME_TAXONOMY_QC.out.gtdbtk_outdir,
                mt_genomes,
                mt_checkm2,
                params.marker_tree_include_references ? REFERENCE_GENOMES.out.gtdbtk_genomes : Channel.value([])
            )
            ch_versions = ch_versions.mix(MARKER_GENE_TREE.out.versions)
        }
    }

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

    if (!params.skip_mobile_elements) {
        MOBILE_ELEMENTS(
            ch_assembly,
            file(params.genomad_db, checkIfExists: true),
            file(params.checkv_db,  checkIfExists: true)
        )
        ch_versions = ch_versions.mix(MOBILE_ELEMENTS.out.versions)
    }

    if (params.run_nonpareil) {
        NONPAREIL(ch_clean)
        ch_versions = ch_versions.mix(NONPAREIL.out.versions)
    }

    // --- Read-stat report (per-sample read tracking across all steps) ---
    READ_STAT_REPORT(
        'metagenome',
        LONG_READ_QC.out.stats.map { meta, stage, t -> t }.collect(),
        ch_scaffold_counts,
        [],
        ch_repmag_abund,
        ch_hq_reps,
        ch_hq_repmag_abund,
        ch_hq_derep_abund,
        ch_hq_ref_abund,
        ch_ws_derep_abund,
        ch_ws_hq_abund
    )
    ch_versions = ch_versions.mix(READ_STAT_REPORT.out.versions)

    // --- Software versions manifest ---
    ch_nfcore_versions = channel.topic('versions')
        .collectFile(name: 'nfcore_versions.yml', newLine: true) { process, tool, version ->
            "\"${process}\":\n    ${tool}: ${version}"
        }
    DUMP_SOFTWARE_VERSIONS(ch_versions.mix(ch_nfcore_versions).collect())
}
