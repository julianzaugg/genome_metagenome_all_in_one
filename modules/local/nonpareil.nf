/*
 * Nonpareil — sequencing depth / coverage redundancy estimation from QC'd reads.
 * Mirrors run_nonpareil.sh (kmer mode on the forward reads).
 */

process NONPAREIL {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.npo"), emit: npo
    path 'versions.yml',                     emit: versions

    script:
    def args  = task.ext.args ?: '-T kmer'
    def fwd   = meta.single_end ? reads : reads[0]
    def ram   = task.memory ? (task.memory.toGiga() * 1000) : 20000
    """
    # nonpareil needs uncompressed fastq input
    seq=${fwd}
    if [[ "\$seq" == *.gz ]]; then zcat \$seq > reads.fastq; seq=reads.fastq; fi

    nonpareil ${args} \\
        -s \$seq \\
        -f fastq \\
        -b ${meta.id} \\
        -R ${ram} \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nonpareil: \$(nonpareil -V 2>&1 | sed 's/Nonpareil v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.npo
    echo '"${task.process}": {nonpareil: stub}' > versions.yml
    """
}
