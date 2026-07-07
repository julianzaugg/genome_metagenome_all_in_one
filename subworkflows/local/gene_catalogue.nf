/*
 * GENE_CATALOGUE — build a dereplicated gene catalogue from per-sample proteins
 * (+ nucleotides) and optionally annotate it with DRAM. Shared by both
 * metagenome tracks.
 *
 *   per-sample .faa + .fna --collect--> CATALOGUE_PREP (complete genes, namespaced)
 *     -> CDHIT (100% + 90%)
 *          -> CATALOGUE_CDS      (nucleotide CDS for 100% reps)
 *          -> CATALOGUE_TABULATE (membership of the 100% catalogue)
 *          -> DRAM (split -> annotate chunks in parallel -> merge; functional
 *                   annotation of the 100% catalogue)
 *
 * When build_expanded is set, an EXPANDED catalogue is ALSO built in the same
 * invocation from the scaffold proteins PLUS external reference-genome proteins,
 * using aliased *_EXPANDED process instances (a process can only be invoked once
 * per workflow graph, so the second catalogue needs distinct names). It gets its
 * own DRAM annotation and is emitted alongside the scaffold-only catalogue.
 */

include { CATALOGUE_PREP; CDHIT; CATALOGUE_CDS; CATALOGUE_TABULATE } from '../../modules/local/gene_catalogue'
include { CATALOGUE_PREP as CATALOGUE_PREP_EXPANDED; CDHIT as CDHIT_EXPANDED;
          CATALOGUE_CDS as CATALOGUE_CDS_EXPANDED; CATALOGUE_TABULATE as CATALOGUE_TABULATE_EXPANDED } from '../../modules/local/gene_catalogue'
include { DRAM_ANNOTATE_SPLIT; DRAM_ANNOTATE; DRAM_ANNOTATE_MERGE } from '../../modules/local/dram'
include { DRAM_ANNOTATE_SPLIT as DRAM_ANNOTATE_SPLIT_EXPANDED; DRAM_ANNOTATE as DRAM_ANNOTATE_EXPANDED;
          DRAM_ANNOTATE_MERGE as DRAM_ANNOTATE_MERGE_EXPANDED } from '../../modules/local/dram'

workflow GENE_CATALOGUE {
    take:
    proteins        // [ meta, faa ] per sample
    nucleotides     // [ meta, fna ] per sample
    identities      // val: comma-separated CD-HIT identities
    dram_db         // path
    run_annotation  // bool
    ref_proteins    // [ meta, faa ] per reference genome (empty channel if none)
    ref_nucleotides // [ meta, fna ] per reference genome (empty channel if none)
    build_expanded  // bool: also build the expanded (scaffold + reference) catalogue

    main:
    ch_faas = proteins.map    { meta, faa -> faa }.collect()
    ch_fnas = nucleotides.map { meta, fna -> fna }.collect()

    CATALOGUE_PREP(ch_faas, ch_fnas)
    CDHIT(CATALOGUE_PREP.out.proteins, identities)
    CATALOGUE_CDS(CDHIT.out.catalogue, CATALOGUE_PREP.out.nucleotides)
    CATALOGUE_TABULATE(CDHIT.out.clusters)

    ch_annotations = Channel.empty()
    ch_versions = Channel.empty()
        .mix(CATALOGUE_PREP.out.versions)
        .mix(CDHIT.out.versions)
        .mix(CATALOGUE_CDS.out.versions)
        .mix(CATALOGUE_TABULATE.out.versions)
    if (run_annotation) {
        // Scatter the catalogue across chunks, annotate in parallel, reassemble.
        DRAM_ANNOTATE_SPLIT(CDHIT.out.catalogue)
        DRAM_ANNOTATE(DRAM_ANNOTATE_SPLIT.out.chunks.flatten(), dram_db)
        DRAM_ANNOTATE_MERGE(
            DRAM_ANNOTATE.out.annotations.collect(),
            DRAM_ANNOTATE.out.annotated_faa.collect()
        )
        ch_annotations = DRAM_ANNOTATE_MERGE.out.annotations
        ch_versions = ch_versions
            .mix(DRAM_ANNOTATE_SPLIT.out.versions)
            .mix(DRAM_ANNOTATE.out.versions.first())
            .mix(DRAM_ANNOTATE_MERGE.out.versions)
    }

    // --- Expanded catalogue: scaffold proteins + reference-genome proteins ---
    ch_exp_catalogue   = Channel.empty()
    ch_exp_cds         = Channel.empty()
    ch_exp_membership  = Channel.empty()
    ch_exp_annotations = Channel.empty()
    if (build_expanded) {
        ch_faas_exp = proteins.mix(ref_proteins).map       { meta, faa -> faa }.collect()
        ch_fnas_exp = nucleotides.mix(ref_nucleotides).map { meta, fna -> fna }.collect()

        CATALOGUE_PREP_EXPANDED(ch_faas_exp, ch_fnas_exp)
        CDHIT_EXPANDED(CATALOGUE_PREP_EXPANDED.out.proteins, identities)
        CATALOGUE_CDS_EXPANDED(CDHIT_EXPANDED.out.catalogue, CATALOGUE_PREP_EXPANDED.out.nucleotides)
        CATALOGUE_TABULATE_EXPANDED(CDHIT_EXPANDED.out.clusters)

        ch_exp_catalogue  = CDHIT_EXPANDED.out.catalogue
        ch_exp_cds        = CATALOGUE_CDS_EXPANDED.out.cds
        ch_exp_membership = CATALOGUE_TABULATE_EXPANDED.out.membership
        ch_versions = ch_versions
            .mix(CATALOGUE_PREP_EXPANDED.out.versions)
            .mix(CDHIT_EXPANDED.out.versions)
            .mix(CATALOGUE_CDS_EXPANDED.out.versions)
            .mix(CATALOGUE_TABULATE_EXPANDED.out.versions)

        if (run_annotation) {
            DRAM_ANNOTATE_SPLIT_EXPANDED(CDHIT_EXPANDED.out.catalogue)
            DRAM_ANNOTATE_EXPANDED(DRAM_ANNOTATE_SPLIT_EXPANDED.out.chunks.flatten(), dram_db)
            DRAM_ANNOTATE_MERGE_EXPANDED(
                DRAM_ANNOTATE_EXPANDED.out.annotations.collect(),
                DRAM_ANNOTATE_EXPANDED.out.annotated_faa.collect()
            )
            ch_exp_annotations = DRAM_ANNOTATE_MERGE_EXPANDED.out.annotations
            ch_versions = ch_versions
                .mix(DRAM_ANNOTATE_SPLIT_EXPANDED.out.versions)
                .mix(DRAM_ANNOTATE_EXPANDED.out.versions.first())
                .mix(DRAM_ANNOTATE_MERGE_EXPANDED.out.versions)
        }
    }

    emit:
    catalogue            = CDHIT.out.catalogue        // primary (first identity)
    catalogues           = CDHIT.out.catalogues       // all identities
    cds                  = CATALOGUE_CDS.out.cds
    membership           = CATALOGUE_TABULATE.out.membership
    annotations          = ch_annotations
    expanded_catalogue   = ch_exp_catalogue
    expanded_cds         = ch_exp_cds
    expanded_membership  = ch_exp_membership
    expanded_annotations = ch_exp_annotations
    versions             = ch_versions
}
