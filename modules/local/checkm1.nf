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
    path 'checkm_qa_results.tsv',         emit: qa
    path 'checkm1',                       emit: outdir
    path 'versions.yml',                  emit: versions

    script:
    def args = task.ext.args ?: ''
    def pplacer_threads = Math.min(task.cpus as int, 40)  // pplacer hangs with many threads (github.com/Ecogenomics/CheckM/issues/341)
    """
    checkm lineage_wf ${args} \\
        -x fasta \\
        -t ${task.cpus} --pplacer_threads ${pplacer_threads} \\
        --tab_table -f checkm_lineage_wf_results.tsv \\
        bins checkm1

    # extended QA table (mirrors run_checkm.sh `checkm qa -o2`)
    checkm qa checkm1/lineage.ms checkm1 \\
        -o 2 --threads ${task.cpus} \\
        --tab_table -f checkm_qa_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$(checkm 2>&1 | grep -m1 -oE 'CheckM v[0-9.]+' | sed 's/CheckM v//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p checkm1
    echo -e "Bin Id\\tCompleteness\\tContamination" > checkm_lineage_wf_results.tsv
    echo -e "Bin Id\\tCompleteness\\tContamination\\tStrain heterogeneity" > checkm_qa_results.tsv
    echo '"${task.process}": {checkm: stub}' > versions.yml
    """
}
