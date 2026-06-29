/*
 * RPKM — SingleM-normalized gene-catalogue abundance from one selected R1 stream.
 */

include {
    RPKM_FILTER_READS;
    RPKM_SINGLEM_MARKERS;
    RPKM_PREBUILT_SINGLEM_MARKERS;
    RPKM_DIAMOND_MAKEDB;
    RPKM_DIAMOND_GENE;
    RPKM_COLLATE_GENE_BLAST;
    RPKM_DIAMOND_SINGLEM;
    RPKM_CALCULATE
} from '../../modules/local/rpkm'

workflow RPKM {
    take:
    reads               // [ meta, selected R1 fastq.gz ] per sample
    catalogue           // gene_catalogue.faa
    singlem_metapackage // SingleM metapackage directory
    singlem_marker_dbs  // optional prebuilt SingleM marker .dmnd directory
    singlem_marker_lengths // optional precomputed SingleM marker lengths TSV
    min_read_length     // val

    main:
    RPKM_FILTER_READS(reads, min_read_length)
    ch_versions = RPKM_FILTER_READS.out.versions
    if (singlem_marker_dbs && singlem_marker_lengths) {
        RPKM_PREBUILT_SINGLEM_MARKERS(singlem_marker_dbs, singlem_marker_lengths)
        ch_marker_dbs = RPKM_PREBUILT_SINGLEM_MARKERS.out.marker_dbs
        ch_marker_lengths = RPKM_PREBUILT_SINGLEM_MARKERS.out.marker_lengths
        ch_versions = ch_versions.mix(RPKM_PREBUILT_SINGLEM_MARKERS.out.versions)
    } else if (!singlem_marker_dbs && !singlem_marker_lengths) {
        RPKM_SINGLEM_MARKERS(singlem_metapackage)
        ch_marker_dbs = RPKM_SINGLEM_MARKERS.out.marker_dbs
        ch_marker_lengths = RPKM_SINGLEM_MARKERS.out.marker_lengths
        ch_versions = ch_versions.mix(RPKM_SINGLEM_MARKERS.out.versions)
    } else {
        error "RPKM requires both --rpkm_singlem_marker_dbs and --rpkm_singlem_marker_lengths when using prebuilt SingleM marker inputs."
    }
    RPKM_DIAMOND_MAKEDB(catalogue)
    ch_versions = ch_versions.mix(RPKM_DIAMOND_MAKEDB.out.versions)

    ch_reads_for_gene = RPKM_FILTER_READS.out.reads
        .combine(RPKM_DIAMOND_MAKEDB.out.db)
        .map { meta, read, db -> [ meta, read, db ] }
    RPKM_DIAMOND_GENE(ch_reads_for_gene)
    RPKM_COLLATE_GENE_BLAST(RPKM_DIAMOND_GENE.out.blast.map { meta, blast -> blast }.collect())
    ch_versions = ch_versions
        .mix(RPKM_DIAMOND_GENE.out.versions)
        .mix(RPKM_COLLATE_GENE_BLAST.out.versions)

    ch_reads_for_singlem = RPKM_FILTER_READS.out.reads
        .combine(ch_marker_dbs)
        .map { meta, read, marker_dbs -> [ meta, read, marker_dbs ] }
    RPKM_DIAMOND_SINGLEM(ch_reads_for_singlem)
    ch_versions = ch_versions.mix(RPKM_DIAMOND_SINGLEM.out.versions)

    RPKM_CALCULATE(
        RPKM_DIAMOND_SINGLEM.out.blast_dir.collect(),
        ch_marker_lengths,
        RPKM_COLLATE_GENE_BLAST.out.blast
    )
    ch_versions = ch_versions.mix(RPKM_CALCULATE.out.versions)

    emit:
    singlem_rpkm    = RPKM_CALCULATE.out.singlem_rpkm
    singlem_means   = RPKM_CALCULATE.out.singlem_means
    normalised_rpkm = RPKM_CALCULATE.out.normalised_rpkm
    mapped_reads    = RPKM_CALCULATE.out.mapped_reads
    gene_blast      = RPKM_COLLATE_GENE_BLAST.out.blast
    marker_blast    = RPKM_DIAMOND_SINGLEM.out.blast_dir
    versions        = ch_versions
}
