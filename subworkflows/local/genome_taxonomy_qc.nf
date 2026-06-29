/*
 * GENOME_TAXONOMY_QC — taxonomy + per-genome analyses on a set of genomes
 * (dereplicated representatives for metagenomes; assemblies for isolates).
 *   GTDB-Tk classify_wf  (all genomes together)
 *   GenomeSPOT           (optional, per-genome: needs proteins)
 *   Barrnap              (optional, per-genome: rRNA / 16S)
 * NOTE: CheckM1/CheckM2 run upstream on the full bin set (they drive dereplication).
 */

include { GTDBTK_CLASSIFYWF }            from '../../modules/nf-core/gtdbtk/classifywf/main'
include { PYRODIGAL as PYRODIGAL_BINS }  from '../../modules/local/pyrodigal'
include { GENOMESPOT; GENOMESPOT_COMBINE } from '../../modules/local/genomespot'
include { BARRNAP }                      from '../../modules/local/barrnap'

workflow GENOME_TAXONOMY_QC {
    take:
    genomes           // collected list of genome fastas
    per_genome        // [ meta, fasta ] per genome
    gtdbtk_db         // path
    genomespot_models // path (may be [])
    run_taxonomy      // bool
    run_genomespot    // bool
    run_barrnap       // bool

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

    ch_genomespot = Channel.empty()
    if (run_genomespot) {
        PYRODIGAL_BINS(per_genome)
        ch_gs_in = per_genome.join(PYRODIGAL_BINS.out.faa)   // [meta, fasta, faa]
        GENOMESPOT(ch_gs_in, genomespot_models)
        ch_genomespot = GENOMESPOT.out.predictions
        GENOMESPOT_COMBINE(GENOMESPOT.out.predictions.map { _meta, tsv -> tsv }.collect())
    }

    ch_rrna = Channel.empty()
    if (run_barrnap) {
        BARRNAP(per_genome)
        ch_rrna = BARRNAP.out.rrna
    }

    emit:
    gtdbtk     = ch_gtdbtk
    genomespot = ch_genomespot
    rrna       = ch_rrna
}
