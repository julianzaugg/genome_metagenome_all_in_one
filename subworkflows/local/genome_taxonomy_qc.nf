/*
 * GENOME_TAXONOMY_QC — taxonomy + per-genome analyses on a set of genomes
 * (all bins for metagenomes; assemblies for isolates).
 *   GTDB-Tk classify_wf  (all genomes together)
 *   GenomeSPOT           (optional, per-genome: needs proteins)
 *   Barrnap              (optional, per-genome: rRNA / 16S)
 * NOTE: CheckM1/CheckM2 run upstream on the full bin set (they drive dereplication).
 */

include { GTDBTK_CLASSIFYWF }            from '../../modules/nf-core/gtdbtk/classifywf/main'
include { PYRODIGAL as PYRODIGAL_BINS }  from '../../modules/local/pyrodigal'
include { GENOMESPOT; GENOMESPOT_COMBINE } from '../../modules/local/genomespot'
include { BARRNAP }                      from '../../modules/local/barrnap'
include { DRAM_ANNOTATE_BINS; DRAM_DISTILL } from '../../modules/local/dram'

workflow GENOME_TAXONOMY_QC {
    take:
    genomes           // collected list of genome fastas
    per_genome        // [ meta, fasta ] per genome
    gtdbtk_db         // path
    genomespot_models // path (may be [])
    dram_db           // path (may be [])
    run_taxonomy      // bool
    run_genomespot    // bool
    run_barrnap       // bool
    run_dram_bins     // bool

    main:
    ch_gtdbtk        = Channel.empty()
    ch_gtdbtk_outdir = Channel.empty()
    ch_versions      = Channel.empty()
    if (run_taxonomy) {
        GTDBTK_CLASSIFYWF(
            genomes.map { g -> [ [id: 'all_genomes'], g ] },
            [ 'gtdb', gtdbtk_db ],
            false
        )
        ch_gtdbtk        = GTDBTK_CLASSIFYWF.out.summary
        ch_gtdbtk_outdir = GTDBTK_CLASSIFYWF.out.gtdb_outdir
        // GTDBTK_CLASSIFYWF emits versions via topic: versions (handled globally)
    }

    // Run pyrodigal once if either GenomeSPOT or DRAM-bins needs per-bin proteins
    ch_proteins = Channel.empty()
    if (run_genomespot || run_dram_bins) {
        PYRODIGAL_BINS(per_genome)
        ch_proteins = PYRODIGAL_BINS.out.faa
        ch_versions = ch_versions.mix(PYRODIGAL_BINS.out.versions)
    }

    ch_genomespot = Channel.empty()
    if (run_genomespot) {
        ch_gs_in = per_genome.join(ch_proteins)   // [meta, fasta, faa]
        GENOMESPOT(ch_gs_in, genomespot_models)
        ch_genomespot = GENOMESPOT.out.predictions
        GENOMESPOT_COMBINE(GENOMESPOT.out.predictions.map { _meta, tsv -> tsv }.collect())
        ch_versions = ch_versions.mix(GENOMESPOT.out.versions)
        ch_versions = ch_versions.mix(GENOMESPOT_COMBINE.out.versions)
    }

    ch_dram_annotations = Channel.empty()
    ch_dram_distilled   = Channel.empty()
    if (run_dram_bins) {
        DRAM_ANNOTATE_BINS(ch_proteins, dram_db)
        DRAM_DISTILL(DRAM_ANNOTATE_BINS.out.annotations, dram_db)
        ch_dram_annotations = DRAM_ANNOTATE_BINS.out.annotations
        ch_dram_distilled   = DRAM_DISTILL.out.distilled
        ch_versions = ch_versions.mix(DRAM_ANNOTATE_BINS.out.versions)
        ch_versions = ch_versions.mix(DRAM_DISTILL.out.versions)
    }

    ch_rrna = Channel.empty()
    if (run_barrnap) {
        BARRNAP(per_genome)
        ch_rrna = BARRNAP.out.rrna
        ch_versions = ch_versions.mix(BARRNAP.out.versions)
    }

    emit:
    gtdbtk           = ch_gtdbtk
    gtdbtk_outdir    = ch_gtdbtk_outdir
    genomespot       = ch_genomespot
    rrna             = ch_rrna
    dram_annotations = ch_dram_annotations
    dram_distilled   = ch_dram_distilled
    versions         = ch_versions
}
