/*
 * LONG_READ_QC — resolve Nanopore FASTQ/POD5 input, adapter trim, length/quality
 * filter, and read stats.
 */

include { DORADO_BASECALL; DORADO_DEMUX; PORECHOP; FASTPLONG } from '../../modules/local/long_reads'
include { FASTQ_GZIP_TEST } from '../../modules/local/validate'
include { SEQKIT_STATS }    from '../../modules/local/read_stats'

workflow LONG_READ_QC {
    take:
    direct_reads        // [ meta, long_reads ] rows used directly
    pod5_for_basecall   // [ meta, pod5_dir ] rows requiring Dorado
    dorado_model        // val
    dorado_barcode_kit  // val
    dorado_device       // val
    run_basecalling     // bool
    run_porechop        // bool
    run_fastplong       // bool

    main:
    ch_raw = direct_reads

    if (run_basecalling) {
        DORADO_BASECALL(pod5_for_basecall, dorado_model, dorado_barcode_kit, dorado_device)
        DORADO_DEMUX(DORADO_BASECALL.out.bam)
        ch_raw = ch_raw.mix(DORADO_DEMUX.out.reads)
    }

    FASTQ_GZIP_TEST(ch_raw)
    ch_reads = FASTQ_GZIP_TEST.out.reads
    ch_read_stats = ch_reads.map { meta, reads -> [ meta, 'raw_long', reads ] }

    if (run_porechop) {
        PORECHOP(ch_reads)
        ch_reads = PORECHOP.out.reads
        ch_read_stats = ch_read_stats.mix(PORECHOP.out.reads.map { meta, reads -> [ meta, 'porechop', reads ] })
    }

    if (run_fastplong) {
        FASTPLONG(ch_reads)
        ch_reads = FASTPLONG.out.reads
        ch_read_stats = ch_read_stats.mix(FASTPLONG.out.reads.map { meta, reads -> [ meta, 'fastplong', reads ] })
    }

    SEQKIT_STATS(ch_read_stats)

    emit:
    reads = ch_reads
    stats = SEQKIT_STATS.out.stats
}
