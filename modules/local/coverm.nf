/*
 * CoverM — bin dereplication (cluster) + read mapping (genome / contig).
 * Mirrors run_bin_dereplication.sh and run_coverm_bins.sh / run_coverm_scaffolds.sh.
 */

process COVERM_CLUSTER {
    label 'process_high'

    input:
    path(bins, stageAs: 'bins/*')
    path(checkm2_report)   // [] if CheckM2 skipped
    path(checkm1_report)   // [] if CheckM1 skipped

    output:
    path 'representatives/*.fasta',                  emit: representatives
    path 'high_quality_representatives/*.fasta',     emit: hq_representatives, optional: true
    path 'cluster_definition.tsv',                   emit: clusters
    path 'versions.yml',                             emit: versions

    script:
    def args    = task.ext.args ?: '--precluster-method finch --ani 0.97'
    def quality = task.ext.quality ?: 50   // completeness - weight*contamination >= quality
    def weight  = task.ext.weight  ?: 3
    """
    # CoverM picks representatives by quality: prefer CheckM2, else CheckM1 tab-table.
    qcflag=""
    if [ -s "${checkm2_report}" ]; then
        qcflag="--checkm2-quality-report ${checkm2_report}"
    elif [ -s "${checkm1_report}" ]; then
        qcflag="--checkm-tab-table ${checkm1_report}"
    fi

    coverm cluster ${args} \\
        -t ${task.cpus} \\
        --genome-fasta-directory bins \\
        --genome-fasta-extension fasta \\
        --output-representative-fasta-directory representatives \\
        --output-cluster-definition cluster_definition.tsv \\
        \$qcflag

    # High-quality subset: a representative is HQ if completeness - ${weight}*contamination
    # >= ${quality} in EITHER CheckM report (union). Column indices found by header name.
    : > hq_ids.txt
    pass_ids() {
        [ -s "\$1" ] || return
        awk -F '\\t' -v q=${quality} -v k=${weight} '
            NR==1 { for(i=1;i<=NF;i++){ if(\$i=="Completeness")cc=i; if(\$i=="Contamination")ct=i;
                                        if(\$i=="Name"||\$i=="Bin Id"||\$i=="genome")id=i } next }
            (cc && ct && id && (\$cc - k*\$ct) >= q) { print \$id }
        ' "\$1" >> hq_ids.txt
    }
    pass_ids "${checkm2_report}"
    pass_ids "${checkm1_report}"

    mkdir -p high_quality_representatives
    sort -u hq_ids.txt | while read -r bin_id; do
        [ -n "\$bin_id" ] && [ -f "representatives/\${bin_id}.fasta" ] && cp "representatives/\${bin_id}.fasta" high_quality_representatives/
    done
    echo "HQ = completeness - ${weight}*contamination >= ${quality} in CheckM1 OR CheckM2" > high_quality_representatives/README.txt

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
