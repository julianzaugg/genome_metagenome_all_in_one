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
    printf "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tN50\tQ20(%%)\tQ30(%%)\tGC(%%)\n" \
        > ${meta.id}.${stage}.seqkit_stats.tsv
    printf "reads.fastq.gz\tFASTQ\tDNA\t1000000\t150000000\t50\t150.0\t151\t150\t150\t151\t150\t95.0\t90.0\t50.0\n" \
        >> ${meta.id}.${stage}.seqkit_stats.tsv
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}
