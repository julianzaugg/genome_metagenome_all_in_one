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
    def args      = task.ext.args ?: '--precluster-method finch --ani 0.95'
    def quality   = task.ext.quality ?: 50   // completeness - weight*contamination >= quality
    def weight    = task.ext.weight  ?: 3
    def hq_source = task.ext.hq_source ?: 'both'   // classify HQ by checkm1, checkm2, or both
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
    # >= ${quality} in the selected CheckM report(s) (${hq_source}). Columns found by header name.
    : > hq_ids.txt
    pass_ids() {
        [ -s "\$1" ] || return
        awk -F '\\t' -v q=${quality} -v k=${weight} '
            NR==1 { for(i=1;i<=NF;i++){ if(\$i=="Completeness")cc=i; if(\$i=="Contamination")ct=i;
                                        if(\$i=="Name"||\$i=="Bin Id"||\$i=="genome")id=i } next }
            (cc && ct && id && (\$cc - k*\$ct) >= q) { print \$id }
        ' "\$1" >> hq_ids.txt
    }
    case "${hq_source}" in
        both)    pass_ids "${checkm2_report}"; pass_ids "${checkm1_report}" ;;
        checkm2) pass_ids "${checkm2_report}" ;;
        checkm1) pass_ids "${checkm1_report}" ;;
    esac

    mkdir -p high_quality_representatives
    sort -u hq_ids.txt | while read -r bin_id; do
        [ -n "\$bin_id" ] && [ -f "representatives/\${bin_id}.fasta" ] && cp "representatives/\${bin_id}.fasta" high_quality_representatives/
    done
    echo "HQ = completeness - ${weight}*contamination >= ${quality} (source: ${hq_source})" > high_quality_representatives/README.txt

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

/*
 * HQ-first dereplication: filter the FULL bin set down to HQ MAGs, THEN cluster.
 * Unlike COVERM_CLUSTER (cluster first, then pick HQ reps), this guarantees every
 * HQ cluster is represented by an HQ genome. representatives = dereplicated HQ set.
 */
process COVERM_CLUSTER_HQ {
    label 'process_high'

    input:
    path(bins, stageAs: 'bins/*')
    path(checkm2_report)   // [] if CheckM2 skipped
    path(checkm1_report)   // [] if CheckM1 skipped

    output:
    path 'representatives/*.fasta', emit: representatives, optional: true
    path 'cluster_definition.tsv',  emit: clusters
    path 'versions.yml',            emit: versions

    script:
    def args      = task.ext.args ?: '--precluster-method finch --ani 0.95'
    def quality   = task.ext.quality ?: 50   // completeness - weight*contamination >= quality
    def weight    = task.ext.weight  ?: 3
    def hq_source = task.ext.hq_source ?: 'both'   // classify HQ by checkm1, checkm2, or both
    """
    # Identify HQ bins first: completeness - ${weight}*contamination >= ${quality} in the
    # selected CheckM report(s) (${hq_source}). Column indices found by header name.
    : > hq_ids.txt
    pass_ids() {
        [ -s "\$1" ] || return
        awk -F '\\t' -v q=${quality} -v k=${weight} '
            NR==1 { for(i=1;i<=NF;i++){ if(\$i=="Completeness")cc=i; if(\$i=="Contamination")ct=i;
                                        if(\$i=="Name"||\$i=="Bin Id"||\$i=="genome")id=i } next }
            (cc && ct && id && (\$cc - k*\$ct) >= q) { print \$id }
        ' "\$1" >> hq_ids.txt
    }
    case "${hq_source}" in
        both)    pass_ids "${checkm2_report}"; pass_ids "${checkm1_report}" ;;
        checkm2) pass_ids "${checkm2_report}" ;;
        checkm1) pass_ids "${checkm1_report}" ;;
    esac

    # Stage only the HQ bins for clustering.
    mkdir -p hq_bins
    sort -u hq_ids.txt | while read -r bin_id; do
        [ -n "\$bin_id" ] && [ -f "bins/\${bin_id}.fasta" ] && cp "bins/\${bin_id}.fasta" hq_bins/
    done
    echo -e "representative\\tmember" > cluster_definition.tsv

    # CoverM picks representatives by quality: prefer CheckM2, else CheckM1 tab-table.
    qcflag=""
    if [ -s "${checkm2_report}" ]; then
        qcflag="--checkm2-quality-report ${checkm2_report}"
    elif [ -s "${checkm1_report}" ]; then
        qcflag="--checkm-tab-table ${checkm1_report}"
    fi

    # No HQ bins passed the filter -> leave representatives empty (emit is optional).
    if ls hq_bins/*.fasta >/dev/null 2>&1; then
        coverm cluster ${args} \\
            -t ${task.cpus} \\
            --genome-fasta-directory hq_bins \\
            --genome-fasta-extension fasta \\
            --output-representative-fasta-directory representatives \\
            --output-cluster-definition cluster_definition.tsv \\
            \$qcflag
    fi

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
    # No reference genomes staged (e.g. no bins passed the HQ filter) -> emit an
    # empty table instead of letting coverm fail on an unmatched fasta glob.
    if ls genomes/*.fasta >/dev/null 2>&1; then
        coverm genome ${args} \\
            --threads ${task.cpus} \\
            --genome-fasta-files genomes/*.fasta \\
            --genome-fasta-extension fasta \\
            --output-file ${meta.id}_abundances.tsv \\
            ${reads_arg}
    else
        echo -e "Genome\\t${meta.id} Relative Abundance (%)" > ${meta.id}_abundances.tsv
    fi

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

process COVERM_CONTIG {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads), path(scaffolds)

    output:
    tuple val(meta), path("${meta.id}_counts.tsv"), emit: counts
    tuple val(meta), path("*.bam"), optional: true, emit: bams
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    def reads_arg = meta.single_end ? "--single ${reads}" : "-c ${reads[0]} ${reads[1]}"
    """
    coverm genome ${args} \\
        --threads ${task.cpus} \\
        --genome-fasta-files ${scaffolds} \\
        --genome-fasta-extension fasta \\
        --output-file ${meta.id}_counts.tsv \\
        ${reads_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/coverm //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "Genome\\t${meta.id} Covered Fraction\\t${meta.id} Mean\\t${meta.id} Count" > ${meta.id}_counts.tsv
    echo '"${task.process}": {coverm: stub}' > versions.yml
    """
}
