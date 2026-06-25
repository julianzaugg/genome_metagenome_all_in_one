/*
 * INPUT_CHECK — parse + validate the samplesheet, emit typed channels.
 *
 * Format/existence are validated against assets/schema_input.json by
 * validateParameters() in main.nf. Here we build meta maps, set the hybrid
 * flags (has_short_reads / has_long_reads / has_pod5), enforce per-mode
 * required columns, and split into downstream channels.
 *
 * Emits (each: [ meta, ... ]):
 *   reads_short : [ meta, [fastq_1, fastq_2] ]   (only rows with short reads)
 *   reads_long  : [ meta, long_reads ]           (only rows with long reads)
 *   pod5        : [ meta, pod5_dir ]             (only rows with a pod5 dir)
 *   host        : [ meta, host_ref ]            (only rows with a host ref)
 *   meta        : [ meta ]                       (all rows; for cross-sample grouping)
 */

def is_set = { v -> v != null && v.toString().trim() && v.toString().trim().toLowerCase() != 'na' }

def build_meta = { row, mode ->
    def meta = [:]
    meta.id    = row.sample.trim()
    meta.group = is_set(row.group) ? row.group.trim() : 'all'

    meta.has_short_reads = is_set(row.fastq_1) && is_set(row.fastq_2)
    meta.has_long_reads  = is_set(row.long_reads)
    meta.has_pod5        = is_set(row.pod5_dir)
    meta.single_end      = false

    def platform = mode.startsWith('illumina') ? 'illumina' : 'nanopore'
    def type     = mode.endsWith('isolate')    ? 'isolate'  : 'metagenome'
    meta.platform = platform
    meta.type     = type

    // --- per-mode required-input checks ---
    if (platform == 'illumina' && !meta.has_short_reads) {
        error "[input_check] Sample '${meta.id}' (mode=${mode}) requires paired fastq_1 + fastq_2."
    }
    if (platform == 'nanopore' && !meta.has_long_reads && !meta.has_pod5) {
        error "[input_check] Sample '${meta.id}' (mode=${mode}) requires long_reads or pod5_dir."
    }
    return meta
}

workflow INPUT_CHECK {
    take:
    samplesheet   // path
    mode          // val: the --mode string

    main:
    ch_rows = Channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row -> [ build_meta(row, mode), row ] }

    reads_short = ch_rows
        .filter { meta, row -> meta.has_short_reads }
        .map    { meta, row -> [ meta, [ file(row.fastq_1, checkIfExists: true),
                                          file(row.fastq_2, checkIfExists: true) ] ] }

    reads_long = ch_rows
        .filter { meta, row -> meta.has_long_reads }
        .map    { meta, row -> [ meta, file(row.long_reads, checkIfExists: true) ] }

    pod5 = ch_rows
        .filter { meta, row -> meta.has_pod5 }
        .map    { meta, row -> [ meta, file(row.pod5_dir, checkIfExists: true) ] }

    host = ch_rows
        .filter { meta, row -> is_set(row.host_ref) }
        .map    { meta, row -> [ meta, file(row.host_ref, checkIfExists: true) ] }

    meta = ch_rows.map { meta, row -> meta }

    emit:
    reads_short
    reads_long
    pod5
    host
    meta
}
