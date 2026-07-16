/*
 * GENOME_TAXONOMY_QC — taxonomy + per-genome analyses on a set of genomes
 * (all bins for metagenomes; assemblies for isolates).
 *   GTDB-Tk classify_wf  (all genomes together)
 *   GenomeSPOT           (optional, per-genome: needs proteins)
 *   Barrnap              (optional, per-genome: rRNA / 16S)
 * NOTE: CheckM1/CheckM2 run upstream on the full bin set (they drive dereplication).
 */

include { GTDBTK_CLASSIFYWF }            from '../../modules/nf-core/gtdbtk/classifywf/main'
include { GTDBTK_RESTORE_NAMES }         from '../../modules/local/gtdbtk_restore'
include { PYRODIGAL as PYRODIGAL_BINS }  from '../../modules/local/pyrodigal'
include { GENOMESPOT; GENOMESPOT_COMBINE } from '../../modules/local/genomespot'
include { BARRNAP }                      from '../../modules/local/barrnap'
include { DRAM_ANNOTATE_BINS; DRAM_COMBINE_ANNOTATIONS; DRAM_DISTILL } from '../../modules/local/dram'

workflow GENOME_TAXONOMY_QC {
    take:
    genomes            // collected list of genome fastas
    per_genome         // [ meta, fasta ] per genome
    gtdbtk_db          // path
    genomespot_models  // path (may be [])
    dram_db            // path (may be [])
    run_taxonomy       // bool
    run_genomespot     // bool
    run_barrnap        // bool
    run_dram_bins      // bool
    reference_genomes  // collected list of USERREF_-prefixed reference fastas (may be [])
    restore_names      // bool: strip the USERREF_ prefix from the published GTDB-Tk output

    main:
    ch_gtdbtk        = Channel.empty()
    ch_gtdbtk_outdir = Channel.empty()
    ch_versions      = Channel.empty()
    if (run_taxonomy) {
        // GTDB-Tk classifies the bins together with the (prefixed) reference genomes
        // in a single run so both are placed in one per-domain tree.
        GTDBTK_CLASSIFYWF(
            genomes.flatten().mix(reference_genomes.flatten()).collect().map { g -> [ [id: 'all_genomes'], g ] },
            [ 'gtdb', gtdbtk_db ],
            false
        )
        ch_gtdbtk_outdir = GTDBTK_CLASSIFYWF.out.gtdb_outdir   // RAW (prefixed) — consumed by the marker tree
        // GTDBTK_CLASSIFYWF emits versions via topic: versions (handled globally)

        if (restore_names) {
            GTDBTK_RESTORE_NAMES(GTDBTK_CLASSIFYWF.out.gtdb_outdir)
            ch_gtdbtk   = GTDBTK_RESTORE_NAMES.out.summary       // published summary with original names
            ch_versions = ch_versions.mix(GTDBTK_RESTORE_NAMES.out.versions)
        } else {
            ch_gtdbtk = GTDBTK_CLASSIFYWF.out.summary
        }
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
        // Concatenate all per-bin annotations (each tagged with its bin via the
        // 'fasta' column) and distill once into a single cross-MAG summary.
        DRAM_COMBINE_ANNOTATIONS(DRAM_ANNOTATE_BINS.out.annotations.map { _meta, tsv -> tsv }.collect())
        DRAM_DISTILL(
            DRAM_COMBINE_ANNOTATIONS.out.annotations.map { tsv -> [ [id: 'all_bins'], tsv ] },
            dram_db
        )
        ch_dram_annotations = DRAM_ANNOTATE_BINS.out.annotations
        ch_dram_distilled   = DRAM_DISTILL.out.distilled
        ch_versions = ch_versions.mix(DRAM_ANNOTATE_BINS.out.versions)
        ch_versions = ch_versions.mix(DRAM_COMBINE_ANNOTATIONS.out.versions)
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
