/*
 * GenomeSPOT — predict growth conditions (temp/pH/salinity/O2) per bin.
 * Mirrors run_genomespot.sh (genome_spot.genome_spot on contigs + proteins).
 */

process GENOMESPOT {
    tag   { meta.id }
    label 'process_low'

    input:
    tuple val(meta), path(contigs), path(proteins)
    path(models)

    output:
    tuple val(meta), path("${meta.id}.predictions.tsv"), emit: predictions
    path 'versions.yml',                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python -m genome_spot.genome_spot ${args} \\
        --contigs ${contigs} \\
        --proteins ${proteins} \\
        --output-prefix ${meta.id} \\
        --models ${models}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomespot: \$(python -m genome_spot.genome_spot --version 2>&1 | tail -1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo -e "target\\tvalue" > ${meta.id}.predictions.tsv
    echo '"${task.process}": {genomespot: stub}' > versions.yml
    """
}
