/*
 * Pyrodigal — gene prediction from assembled scaffolds.
 * Mirrors run_pyrodigal_scaffolds.sh (meta mode; emits .faa/.fna/.gff).
 */

process PYRODIGAL {
    tag   { meta.id }
    label 'process_low'

    input:
    tuple val(meta), path(scaffolds)

    output:
    tuple val(meta), path("${meta.id}.faa"), emit: faa
    tuple val(meta), path("${meta.id}.fna"), emit: fna
    tuple val(meta), path("${meta.id}.gff"), emit: gff
    path 'versions.yml',                     emit: versions

    script:
    def args = task.ext.args ?: '-p meta'
    """
    pyrodigal -i ${scaffolds} ${args} \\
        -a ${meta.id}.faa \\
        -d ${meta.id}.fna \\
        -f gff > ${meta.id}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyrodigal: \$(pyrodigal --version 2>&1 | sed 's/pyrodigal //')
    END_VERSIONS
    """

    stub:
    """
    echo ">${meta.id}_1 # partial=00" > ${meta.id}.faa; echo "MAAA" >> ${meta.id}.faa
    echo ">${meta.id}_1" > ${meta.id}.fna; echo "ATGGCT" >> ${meta.id}.fna
    echo "##gff-version 3" > ${meta.id}.gff
    echo '"${task.process}": {pyrodigal: stub}' > versions.yml
    """
}
