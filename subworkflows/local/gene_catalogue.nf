/*
 * GENE_CATALOGUE — build a dereplicated gene catalogue from per-sample proteins
 * (+ nucleotides) and optionally annotate it with DRAM. Shared by both
 * metagenome tracks.
 *
 *   per-sample .faa + .fna --collect--> CATALOGUE_PREP (complete genes, namespaced)
 *     -> CDHIT (100% + 90%)
 *          -> CATALOGUE_CDS      (nucleotide CDS for 100% reps)
 *          -> CATALOGUE_TABULATE (membership of the 100% catalogue)
 *          -> DRAM_ANNOTATE      (functional annotation of the 100% catalogue)
 */

include { CATALOGUE_PREP; CDHIT; CATALOGUE_CDS; CATALOGUE_TABULATE } from '../../modules/local/gene_catalogue'
include { DRAM_ANNOTATE }                                          from '../../modules/local/dram'

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
    if (run_annotation) {
        DRAM_ANNOTATE(CDHIT.out.catalogue, dram_db)
        ch_annotations = DRAM_ANNOTATE.out.annotations
    }

    emit:
    catalogue   = CDHIT.out.catalogue        // primary (first identity)
    catalogues  = CDHIT.out.catalogues       // all identities
    cds         = CATALOGUE_CDS.out.cds
    membership  = CATALOGUE_TABULATE.out.membership
    annotations = ch_annotations
}
