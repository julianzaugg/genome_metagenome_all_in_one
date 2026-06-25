/*
 * Host / contaminant read removal.
 *  - MINIMAP2_HOSTFILTER : paired short reads (mirrors run_filter_human.sh:
 *        minimap2 -ax sr | samtools fastq -f 13  → keep read pairs where neither
 *        mate mapped to the host).
 *  - CLEANIFIER          : long / single reads (mirrors run_cleanifier.sh).
 * The HOST_REMOVAL subworkflow picks one per-sample based on meta.single_end.
 */

process MINIMAP2_HOSTFILTER {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(host_ref)

    output:
    tuple val(meta), path("${meta.id}.hostfilt_{1,2}.fastq.gz"), emit: reads
    path 'versions.yml',                                         emit: versions

    script:
    def args = task.ext.args ?: '-ax sr'
    """
    minimap2 ${args} -t ${task.cpus} ${host_ref} ${reads[0]} ${reads[1]} \\
        | samtools fastq --threads ${task.cpus} -n -f 13 - \\
            -1 ${meta.id}.hostfilt_1.fastq.gz \\
            -2 ${meta.id}.hostfilt_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.hostfilt_1.fastq.gz
    echo | gzip > ${meta.id}.hostfilt_2.fastq.gz
    echo '"${task.process}": {minimap2: stub}' > versions.yml
    """
}

process CLEANIFIER {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads
    path 'versions.yml',                          emit: versions

    script:
    def args = task.ext.args ?: '--compression gz --compression-level 7 --threshold 0.5 --buffersize 20'
    """
    cleanifier filter ${args} \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        --index ${index} \\
        --prefix ${meta.id}
    mv ${meta.id}_keep.fastq.gz ${meta.id}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cleanifier: \$(cleanifier --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.fastq.gz
    echo '"${task.process}": {cleanifier: stub}' > versions.yml
    """
}
