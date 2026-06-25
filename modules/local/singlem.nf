/*
 * SingleM — marker-gene taxonomic profiling from raw reads.
 * Mirrors run_singlem.sh (singlem pipe). Run per-sample; merge downstream.
 */

process SINGLEM_PIPE {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(metapackage)

    output:
    tuple val(meta), path("${meta.id}.condensed.tsv"),       emit: profile
    tuple val(meta), path("${meta.id}.otu_table.tsv"),       emit: otu
    tuple val(meta), path("${meta.id}.otu_table.archive"),   emit: archive
    path 'versions.yml',                                     emit: versions

    script:
    def args  = task.ext.args ?: ''
    def input = meta.single_end ? "--forward ${reads}" : "--forward ${reads[0]} --reverse ${reads[1]}"
    """
    singlem pipe ${args} \\
        --metapackage ${metapackage} \\
        ${input} \\
        --archive-otu-table ${meta.id}.otu_table.archive \\
        --taxonomic-profile ${meta.id}.condensed.tsv \\
        --otu-table ${meta.id}.otu_table.tsv \\
        --threads ${task.cpus} \\
        --assignment-threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.condensed.tsv ${meta.id}.otu_table.tsv ${meta.id}.otu_table.archive
    echo '"${task.process}": {singlem: stub}' > versions.yml
    """
}
