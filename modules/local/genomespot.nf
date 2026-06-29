/*
 * GenomeSPOT — predict growth conditions (temp/pH/salinity/O2) per bin.
 * Mirrors run_genomespot.sh (genome_spot.genome_spot on contigs + proteins).
 */

process GENOMESPOT_COMBINE {
    label 'process_single'

    input:
    path(tsvs, stageAs: 'predictions/*')

    output:
    path 'genomespot_combined.tsv', emit: combined
    path 'versions.yml',            emit: versions

    script:
    """
    OUTPUT=genomespot_combined.tsv
    header_written=0
    for f in predictions/*.tsv; do
        bin_id=\$(basename "\$f" .predictions.tsv)
        if [[ \$header_written -eq 0 ]]; then
            awk -v id="\$bin_id" 'NR==1{print "bin_id\\t"\$0; next}{print id"\\t"\$0}' "\$f" > "\$OUTPUT"
            header_written=1
        else
            awk -v id="\$bin_id" 'NR>1{print id"\\t"\$0}' "\$f" >> "\$OUTPUT"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomespot: \$(python -m genome_spot.genome_spot --version 2>&1 | tail -1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo -e "bin_id\\ttarget\\tvalue" > genomespot_combined.tsv
    echo '"${task.process}": {genomespot: stub}' > versions.yml
    """
}

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
