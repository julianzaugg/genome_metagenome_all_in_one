/*
 * Host / contaminant read removal with cleanifier (k-mer based).
 * Mirrors run_cleanifier.sh; paired-end uses --fastq R1 --pairs R2 per
 * https://gitlab.com/rahmannlab/cleanifier . Kept (host-removed) reads are
 * renamed to standard names for downstream steps.
 */

process CLEANIFIER {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta.id}.clean*.fastq.gz"), emit: reads
    path 'versions.yml',                                 emit: versions

    script:
    def args  = task.ext.args ?: '--compression gz --compression-level 7 --threshold 0.5 --buffersize 20'
    def input = meta.single_end ? "--fastq ${reads}" : "--fastq ${reads[0]} --pairs ${reads[1]}"
    """
    cleanifier filter ${args} \\
        --threads ${task.cpus} \\
        --index ${index} \\
        ${input} \\
        --prefix ${meta.id}

    # cleanifier writes kept reads as <prefix>...keep...{,.1/.2}.f(ast)q.gz —
    # standardise names (adjust the glob if your cleanifier version differs).
    if [ "${meta.single_end}" = "true" ]; then
        mv \$(ls ${meta.id}*keep*.f*q.gz | head -1) ${meta.id}.clean.fastq.gz
    else
        keeps=(\$(ls ${meta.id}*keep*.f*q.gz | sort))
        mv "\${keeps[0]}" ${meta.id}.clean_1.fastq.gz
        mv "\${keeps[1]}" ${meta.id}.clean_2.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cleanifier: \$(cleanifier --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    if [ "${meta.single_end}" = "true" ]; then
        echo | gzip > ${meta.id}.clean.fastq.gz
    else
        echo | gzip > ${meta.id}.clean_1.fastq.gz
        echo | gzip > ${meta.id}.clean_2.fastq.gz
    fi
    echo '"${task.process}": {cleanifier: stub}' > versions.yml
    """
}
