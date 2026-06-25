/*
 * CheckM1 — lineage_wf bin completeness/contamination (run alongside CheckM2).
 * Mirrors run_checkm.sh. Operates on the full set of bins at once.
 */

process CHECKM1_LINEAGEWF {
    label 'process_high'

    input:
    path(bins, stageAs: 'bins/*')

    output:
    path 'checkm_lineage_wf_results.tsv', emit: summary
    path 'checkm1',                       emit: outdir
    path 'versions.yml',                  emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    checkm lineage_wf ${args} \\
        -x fasta \\
        -t ${task.cpus} \\
        --tab_table -f checkm_lineage_wf_results.tsv \\
        bins checkm1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$(checkm 2>&1 | grep -m1 -oE 'CheckM v[0-9.]+' | sed 's/CheckM v//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p checkm1
    echo -e "Bin Id\\tCompleteness\\tContamination" > checkm_lineage_wf_results.tsv
    echo '"${task.process}": {checkm: stub}' > versions.yml
    """
}
