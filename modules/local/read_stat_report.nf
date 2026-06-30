/*
 * Aggregate per-sample read counts from every processing step into one TSV.
 *
 * Sources are the pipeline's existing per-sample outputs: SeqKit stats per QC
 * stage and CoverM read-mapping tables. Inputs are staged into distinct
 * directories (the ONT and short-read CoverM contig tables share the
 * `{id}_counts.tsv` name) and each optional input defaults to [].
 *   report_mode: 'metagenome' | 'isolate' — selects the column layout.
 */

process READ_STAT_REPORT {
    label 'process_single'

    input:
    val(report_mode)
    path(seqkit_stats,   stageAs: 'seqkit/*')
    path(scaffold_counts, stageAs: 'scaffolds/*')
    path(scaffold_sr,    stageAs: 'scaffolds_sr/*')
    path(repmag_abund,   stageAs: 'repmags/*')
    path(hq_reps,        stageAs: 'hq/*')

    output:
    path 'read_stat_report.tsv', emit: report
    path 'versions.yml',         emit: versions

    script:
    """
    build_read_stat_report.py \\
        --mode ${report_mode} \\
        --seqkit-dir seqkit \\
        --scaffold-dir scaffolds \\
        --scaffold-sr-dir scaffolds_sr \\
        --repmag-dir repmags \\
        --hq-dir hq \\
        --out read_stat_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    build_read_stat_report.py \\
        --mode ${report_mode} \\
        --seqkit-dir seqkit \\
        --scaffold-dir scaffolds \\
        --scaffold-sr-dir scaffolds_sr \\
        --repmag-dir repmags \\
        --hq-dir hq \\
        --out read_stat_report.tsv

    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
