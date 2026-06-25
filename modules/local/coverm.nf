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
    path 'representatives/*.fasta',     emit: representatives
    path 'cluster_definition.tsv',      emit: clusters
    path 'versions.yml',                emit: versions

    script:
    def args = task.ext.args ?: '--precluster-method finch --ani 0.97'
    def qc   = checkm2_report ? "--checkm2-quality-report ${checkm2_report}" : ''
    """
    coverm cluster ${args} \\
        -t ${task.cpus} \\
        --genome-fasta-directory bins \\
        --genome-fasta-extension fasta \\
        --output-representative-fasta-directory representatives \\
        --output-cluster-definition cluster_definition.tsv \\
        ${qc}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/coverm //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p representatives
    cp bins/* representatives/ 2>/dev/null || (echo ">c" > representatives/rep.1.fasta; echo "ACGT" >> representatives/rep.1.fasta)
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
