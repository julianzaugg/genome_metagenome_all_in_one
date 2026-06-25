/*
 * GENOME_TAXONOMY_QC — taxonomy + quality on a set of genomes (dereplicated bins
 * for metagenomes, assemblies for isolates).
 *   GTDB-Tk classify_wf  (all genomes together)
 *   CheckM1 lineage_wf    (optional, run alongside CheckM2 which runs upstream)
 *   GenomeSPOT            (optional, per-genome: needs proteins)
 */

include { GTDBTK_CLASSIFYWF }            from '../../modules/nf-core/gtdbtk/classifywf/main'
include { CHECKM1_LINEAGEWF }            from '../../modules/local/checkm1'
include { PYRODIGAL as PYRODIGAL_BINS }  from '../../modules/local/pyrodigal'
include { GENOMESPOT }                   from '../../modules/local/genomespot'

workflow GENOME_TAXONOMY_QC {
    take:
    genomes          // collected list of genome fastas
    per_genome       // [ meta, fasta ] per genome (for genomespot)
    gtdbtk_db        // path
    genomespot_models // path (may be [])
    run_taxonomy     // bool
    run_checkm1      // bool
    run_genomespot   // bool

    main:
    ch_gtdbtk = Channel.empty()
    if (run_taxonomy) {
        GTDBTK_CLASSIFYWF(
            genomes.map { g -> [ [id: 'all_genomes'], g ] },
            [ 'gtdb', gtdbtk_db ],
            false
        )
        ch_gtdbtk = GTDBTK_CLASSIFYWF.out.summary
    }

    if (run_checkm1) {
        CHECKM1_LINEAGEWF(genomes)
    }

    ch_genomespot = Channel.empty()
    if (run_genomespot) {
        PYRODIGAL_BINS(per_genome)
        ch_gs_in = per_genome.join(PYRODIGAL_BINS.out.faa)   // [meta, fasta, faa]
        GENOMESPOT(ch_gs_in, genomespot_models)
        ch_genomespot = GENOMESPOT.out.predictions
    }

    emit:
    gtdbtk    = ch_gtdbtk
    genomespot = ch_genomespot
}
