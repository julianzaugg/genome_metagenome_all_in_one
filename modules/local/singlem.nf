/*
 * SingleM — marker-gene taxonomic profiling from raw reads.
 * Mirrors run_singlem.sh: one multi-sample singlem pipe call over all samples.
 */

process SINGLEM_PIPE {
    label 'process_high'

    input:
    path(forward_reads)
    path(reverse_reads)
    path(metapackage)

    output:
    path 'metagenome.condensed.tsv',             emit: profile
    path 'metagenome.otu_table.tsv',             emit: otu
    path 'metagenome.otu_table.archive',         emit: archive
    path 'metagenome.prokaryotic_fraction.tsv',  emit: prokaryotic_fraction
    path 'versions.yml',                         emit: versions

    script:
    def args = task.ext.args ?: ''
    def fwd = forward_reads instanceof Collection ? forward_reads : [forward_reads]
    def rev = reverse_reads instanceof Collection ? reverse_reads.findAll { it } : (reverse_reads ? [reverse_reads] : [])
    def input = rev ? "--forward ${fwd.join(' ')} --reverse ${rev.join(' ')}" : "--forward ${fwd.join(' ')}"
    """
    singlem pipe ${args} \\
        --metapackage ${metapackage} \\
        ${input} \\
        --archive-otu-table metagenome.otu_table.archive \\
        --taxonomic-profile metagenome.condensed.tsv \\
        --otu-table metagenome.otu_table.tsv \\
        --threads ${task.cpus} \\
        --assignment-threads ${task.cpus}

    # prokaryotic (microbial) fraction — SingleM >=0.21 renamed microbial_fraction
    # to prokaryotic_fraction; passing the reads makes it COMPUTE the metagenome
    # size (no precomputed --input-metagenome-sizes needed).
    # https://wwood.github.io/singlem/tools/prokaryotic_fraction
    singlem prokaryotic_fraction \\
        --input-profile metagenome.condensed.tsv \\
        --metapackage ${metapackage} \\
        ${input} \\
        --output-tsv metagenome.prokaryotic_fraction.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        singlem: \$(singlem --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch metagenome.condensed.tsv metagenome.otu_table.tsv metagenome.otu_table.archive
    touch metagenome.prokaryotic_fraction.tsv
    echo '"${task.process}": {singlem: stub}' > versions.yml
    """
}
