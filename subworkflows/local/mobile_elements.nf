/*
 * MOBILE_ELEMENTS — virus/plasmid prediction, QC and ANI clustering.
 *   geNomad end_to_end (per sample assembly)
 *     -> pool virus/plasmid sequences across samples (GENOMAD_PROCESS)
 *     -> CheckV end_to_end on pooled viruses -> collate viruses+proviruses
 *     -> ANI clustering (blast -> anicalc -> aniclust) for virus AND plasmid.
 */

include { GENOMAD_ENDTOEND }                from '../../modules/nf-core/genomad/endtoend/main'
include { GENOMAD_PROCESS }                 from '../../modules/local/genomad_process'
include { CHECKV_ENDTOEND }                 from '../../modules/nf-core/checkv/endtoend/main'
include { CHECKV_COLLATE; CHECKV_CLUSTER }  from '../../modules/local/checkv_cluster'

workflow MOBILE_ELEMENTS {
    take:
    assemblies   // [ meta, fasta ] per sample
    genomad_db   // path
    checkv_db    // path

    main:
    GENOMAD_ENDTOEND(assemblies, genomad_db)

    // collect each geNomad per-sample output across samples (drop meta)
    def justfiles = { ch -> ch.map { meta, f -> f }.collect() }
    GENOMAD_PROCESS(
        justfiles(GENOMAD_ENDTOEND.out.virus_fasta),
        justfiles(GENOMAD_ENDTOEND.out.plasmid_fasta),
        justfiles(GENOMAD_ENDTOEND.out.virus_proteins),
        justfiles(GENOMAD_ENDTOEND.out.plasmid_proteins),
        justfiles(GENOMAD_ENDTOEND.out.virus_genes),
        justfiles(GENOMAD_ENDTOEND.out.plasmid_genes),
        justfiles(GENOMAD_ENDTOEND.out.virus_summary),
        justfiles(GENOMAD_ENDTOEND.out.plasmid_summary)
    )

    // CheckV on pooled viral sequences, then collate viruses + proviruses.
    CHECKV_ENDTOEND(GENOMAD_PROCESS.out.virus.map { f -> [ [id: 'pooled_virus'], f ] }, checkv_db)
    ch_collate_in = CHECKV_ENDTOEND.out.viruses.join(CHECKV_ENDTOEND.out.proviruses)
    CHECKV_COLLATE(ch_collate_in)

    // Cluster virus (post-CheckV) and plasmid (geNomad) sequences separately.
    ch_cluster_in = CHECKV_COLLATE.out.fasta.map { f -> ['virus', f] }
        .mix( GENOMAD_PROCESS.out.plasmid.map { f -> ['plasmid', f] } )
    CHECKV_CLUSTER(ch_cluster_in)

    emit:
    clusters = CHECKV_CLUSTER.out.clusters
}
