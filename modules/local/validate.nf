/*
 * Lightweight input validation helpers.
 */

process FASTQ_GZIP_TEST {
    tag   { meta.id }
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), emit: reads
    path 'versions.yml',          emit: versions

    script:
    def files = reads instanceof Collection ? reads : [reads]
    def quoted_files = files.collect { "\"${it}\"" }.join(' ')
    def file_names = files.collect { it.name }.join(' ')
    """
    set -euo pipefail

    if ! gzip -t ${quoted_files}; then
        echo "[FASTQ_GZIP_TEST] Sample ${meta.id} has a corrupt or truncated gzip FASTQ: ${file_names}" >&2
        echo "[FASTQ_GZIP_TEST] Re-copy/re-download the file(s), then resume the pipeline." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: stub
    END_VERSIONS
    """
}
