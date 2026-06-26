/*
 * CoverM — bin dereplication (cluster) + read mapping (genome / contig).
 * Mirrors run_bin_dereplication.sh and run_coverm_bins.sh / run_coverm_scaffolds.sh.
 */

process COVERM_CLUSTER {
    label 'process_high'

    input:
    path(bins, stageAs: 'bins/*')
    path(checkm2_report)

    output:
    path 'representatives/*.fasta',                  emit: representatives
    path 'high_quality_representatives/*.fasta',     emit: hq_representatives, optional: true
    path 'cluster_definition.tsv',                   emit: clusters
    path 'versions.yml',                             emit: versions

    script:
    def args    = task.ext.args ?: '--precluster-method finch --ani 0.97'
    def quality = task.ext.quality ?: 50   // completeness - 3*contamination >= quality
    def qc      = checkm2_report ? "--checkm2-quality-report ${checkm2_report}" : ''
    """
    coverm cluster ${args} \\
        -t ${task.cpus} \\
        --genome-fasta-directory bins \\
        --genome-fasta-extension fasta \\
        --output-representative-fasta-directory representatives \\
        --output-cluster-definition cluster_definition.tsv \\
        ${qc}

    # high-quality subset: completeness - 3*contamination >= ${quality}
    # (mirrors run_bin_dereplication.sh; CheckM2 columns: 1=Name 2=Completeness 3=Contamination)
    mkdir -p high_quality_representatives
    if [ -n "${qc}" ]; then
        awk -F '\\t' 'NR>1 && (\$2 - (\$3*3)) >= ${quality} {print \$1}' ${checkm2_report} \\
        | while read -r bin_id; do
            [ -f "representatives/\${bin_id}.fasta" ] && cp "representatives/\${bin_id}.fasta" high_quality_representatives/
        done
        echo "Representatives with completeness - 3*contamination >= ${quality}" > high_quality_representatives/README.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/coverm //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p representatives high_quality_representatives
    cp bins/* representatives/ 2>/dev/null || (echo ">c" > representatives/rep.1.fasta; echo "ACGT" >> representatives/rep.1.fasta)
    cp representatives/*.fasta high_quality_representatives/ 2>/dev/null || true
    echo -e "representative\\tmember" > cluster_definition.tsv
    echo '"${task.process}": {coverm: stub}' > versions.yml
    """
}

process COVERM_GENOME {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(genomes, stageAs: 'genomes/*')

    output:
    tuple val(meta), path("${meta.id}_abundances.tsv"), emit: abundance
    path 'versions.yml',                                emit: versions

    script:
    def args = task.ext.args ?: ''
    def reads_arg = meta.single_end ? "--single ${reads}" : "-c ${reads[0]} ${reads[1]}"
    """
    coverm genome ${args} \\
        --threads ${task.cpus} \\
        --genome-fasta-files genomes/*.fasta \\
        --genome-fasta-extension fasta \\
        --output-file ${meta.id}_abundances.tsv \\
        ${reads_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/coverm //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "Genome\\t${meta.id} Relative Abundance (%)" > ${meta.id}_abundances.tsv
    echo '"${task.process}": {coverm: stub}' > versions.yml
    """
}
