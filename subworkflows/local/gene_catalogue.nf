/*
 * GENE_CATALOGUE — build a dereplicated gene catalogue from per-sample proteins
 * and (optionally) annotate it with DRAM. Shared by both metagenome tracks.
 *
 *   per-sample .faa  --collect-->  CATALOGUE_PREP (complete genes, namespaced)
 *     -> CDHIT (cluster)  -> CATALOGUE_TABULATE (membership)
 *                         -> DRAM_ANNOTATE (functional annotation)
 */

include { CATALOGUE_PREP; CDHIT; CATALOGUE_TABULATE } from '../../modules/local/gene_catalogue'
include { DRAM_ANNOTATE }                            from '../../modules/local/dram'

workflow GENE_CATALOGUE {
    take:
    proteins      // [ meta, faa ] per sample
    dram_db       // path
    run_annotation // bool

    main:
    ch_faas = proteins.map { meta, faa -> faa }.collect()

    CATALOGUE_PREP(ch_faas)
    CDHIT(CATALOGUE_PREP.out.proteins)
    CATALOGUE_TABULATE(CDHIT.out.clusters)

    ch_annotations = Channel.empty()
    if (run_annotation) {
        DRAM_ANNOTATE(CDHIT.out.catalogue, dram_db)
        ch_annotations = DRAM_ANNOTATE.out.annotations
    }

    emit:
    catalogue   = CDHIT.out.catalogue
    membership  = CATALOGUE_TABULATE.out.membership
    annotations = ch_annotations
}
