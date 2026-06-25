/*
 * ANI-based clustering of viral / plasmid sequences (NOT MCL).
 * Mirrors run_dereplicate_checkv.sh: makeblastdb + all-vs-all blastn -> anicalc.py
 * -> aniclust.py at MIUVIG thresholds (95% ANI / 85% AF). `label` distinguishes
 * the virus vs plasmid run.
 */

process CHECKV_COLLATE {
    // Combine CheckV viruses.fna + proviruses.fna into one fasta for clustering
    // (mirrors `cat viruses.fna proviruses.fna > checkv_viruses.fasta`).
    label 'process_single'

    input:
    tuple val(meta), path(viruses), path(proviruses)

    output:
    path 'checkv_viruses.fasta', emit: fasta

    script:
    """
    cat ${viruses} ${proviruses} > checkv_viruses.fasta
    """

    stub:
    """
    cat ${viruses} ${proviruses} > checkv_viruses.fasta 2>/dev/null || echo ">v" > checkv_viruses.fasta
    """
}

process CHECKV_CLUSTER {
    tag   { label }
    label 'process_high'

    input:
    tuple val(label), path(fasta)

    output:
    tuple val(label), path("${label}_clusters.tsv"), emit: clusters
    tuple val(label), path("${label}_ani.tsv"),      emit: ani
    path 'versions.yml',                             emit: versions

    script:
    def args = task.ext.args ?: '--min_ani 95 --min_tcov 85 --min_qcov 0'
    """
    makeblastdb -in ${fasta} -dbtype nucl -out ${label}_db
    blastn -query ${fasta} -db ${label}_db \\
        -outfmt '6 std qlen slen' -max_target_seqs 10000 \\
        -num_threads ${task.cpus} -out ${label}_blast.tsv

    anicalc.py  -i ${label}_blast.tsv -o ${label}_ani.tsv
    aniclust.py --fna ${fasta} --ani ${label}_ani.tsv --out ${label}_clusters.tsv ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -1 | sed 's/blastn: //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "representative\\tmembers" > ${label}_clusters.tsv
    echo -e "qname\\ttname\\tnum_alns\\tani\\tqcov\\ttcov" > ${label}_ani.tsv
    echo '"${task.process}": {blast: stub}' > versions.yml
    """
}
