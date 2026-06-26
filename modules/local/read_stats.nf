/*
 * SeqKit read statistics for raw and QC-step FASTQ outputs.
 */

process SEQKIT_STATS {
    tag   { "${meta.id}:${stage}" }
    label 'process_single'

    input:
    tuple val(meta), val(stage), path(reads)

    output:
    tuple val(meta), val(stage), path("${meta.id}.${stage}.seqkit_stats.tsv"), emit: stats
    path 'versions.yml', emit: versions

    script:
    def files = reads instanceof Collection ? reads : [reads]
    def quoted_files = files.collect { "\"${it}\"" }.join(' ')
    """
    seqkit stats --tabular --all ${quoted_files} > ${meta.id}.${stage}.seqkit_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/^seqkit //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len" > ${meta.id}.${stage}.seqkit_stats.tsv
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}
