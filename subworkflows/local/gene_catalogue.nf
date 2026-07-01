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
 */

include { CATALOGUE_PREP; CDHIT; CATALOGUE_CDS; CATALOGUE_TABULATE } from '../../modules/local/gene_catalogue'
include { DRAM_ANNOTATE_SPLIT; DRAM_ANNOTATE; DRAM_ANNOTATE_MERGE } from '../../modules/local/dram'

workflow GENE_CATALOGUE {
    take:
    proteins       // [ meta, faa ] per sample
    nucleotides    // [ meta, fna ] per sample
    identities     // val: comma-separated CD-HIT identities
    dram_db        // path
    run_annotation // bool

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

    emit:
    catalogue   = CDHIT.out.catalogue        // primary (first identity)
    catalogues  = CDHIT.out.catalogues       // all identities
    cds         = CATALOGUE_CDS.out.cds
    membership  = CATALOGUE_TABULATE.out.membership
    annotations = ch_annotations
    versions    = ch_versions
}
