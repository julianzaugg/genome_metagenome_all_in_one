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

/*
 * HQ-first dereplication WITH external reference genomes.
 * Like COVERM_CLUSTER_HQ (filter the full bin set to HQ MAGs, then cluster), but the
 * supplied reference genomes are added to the clustering pool UNCONDITIONALLY (they
 * bypass the HQ filter — always included). Representative selection is quality-driven
 * across MAGs and refs via a combined CheckM2 report (MAG report + reference report).
 * representatives = the dereplicated (HQ MAGs + references) set.
 */
process COVERM_CLUSTER_HQ_REF {
    label 'process_high'

    input:
    path(bins, stageAs: 'bins/*')
    path(refs, stageAs: 'refs/*')
    path(checkm2_report)      // [] if CheckM2 skipped
    path(checkm1_report)      // [] if CheckM1 skipped
    path(ref_checkm2_report)  // CheckM2 report for the reference genomes

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

    # Stage the HQ bins for clustering.
    mkdir -p hq_bins
    sort -u hq_ids.txt | while read -r bin_id; do
        [ -n "\$bin_id" ] && [ -f "bins/\${bin_id}.fasta" ] && cp "bins/\${bin_id}.fasta" hq_bins/
    done

    # Add all reference genomes unconditionally. Fail loudly on a name clash with a bin
    # (would make both the staged FASTA and the combined CheckM2 report ambiguous).
    for r in refs/*.fasta; do
        [ -e "\$r" ] || continue
        rb=\$(basename "\$r")
        if [ -f "bins/\$rb" ]; then
            echo "[COVERM_CLUSTER_HQ_REF] Reference genome '\$rb' collides with a bin of the same name. Rename the reference genome." >&2
            exit 1
        fi
        cp "\$r" hq_bins/
    done
    echo -e "representative\\tmember" > cluster_definition.tsv

    # Combined CheckM2 report for representative selection = MAG report + reference report
    # (single header). References were scored with CheckM2 upstream.
    : > combined_checkm2.tsv
    if [ -s "${checkm2_report}" ]; then cat "${checkm2_report}" > combined_checkm2.tsv; fi
    if [ -s "${ref_checkm2_report}" ]; then
        if [ -s combined_checkm2.tsv ]; then
            tail -n +2 "${ref_checkm2_report}" >> combined_checkm2.tsv
        else
            cat "${ref_checkm2_report}" > combined_checkm2.tsv
        fi
    fi

    qcflag=""
    if [ -s combined_checkm2.tsv ]; then
        qcflag="--checkm2-quality-report combined_checkm2.tsv"
    elif [ -s "${checkm1_report}" ]; then
        qcflag="--checkm-tab-table ${checkm1_report}"
    fi

    # No genomes staged (no HQ bins and no refs) -> leave representatives empty (emit optional).
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
    cp bins/* representatives/ 2>/dev/null || true
    cp refs/* representatives/ 2>/dev/null || true
    ls representatives/*.fasta >/dev/null 2>&1 || (echo ">c" > representatives/rep.1.fasta; echo "ACGT" >> representatives/rep.1.fasta)
    echo -e "representative\\tmember" > cluster_definition.tsv
    echo '"${task.process}": {coverm: stub}' > versions.yml
    """
}

/*
 * Within-sample (or within-group) dereplication — a per-unit copy of COVERM_CLUSTER.
 * Identical logic, but carries meta so each sample/group's bins are clustered on
 * their own and the per-unit identity survives for publishing + read mapping.
 */
process COVERM_CLUSTER_WS {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(bins, stageAs: 'bins/*')
    path(checkm2_report)   // [] if CheckM2 skipped
    path(checkm1_report)   // [] if CheckM1 skipped

    output:
    tuple val(meta), path('representatives/*.fasta'),              emit: representatives
    tuple val(meta), path('high_quality_representatives/*.fasta'), emit: hq_representatives, optional: true
    tuple val(meta), path('cluster_definition.tsv'),              emit: clusters
    path 'versions.yml',                                          emit: versions

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
 * Within-sample (or within-group) HQ-first dereplication — a per-unit copy of
 * COVERM_CLUSTER_HQ: filter the unit's bins down to HQ MAGs, THEN cluster.
 */
process COVERM_CLUSTER_HQ_WS {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(bins, stageAs: 'bins/*')
    path(checkm2_report)   // [] if CheckM2 skipped
    path(checkm1_report)   // [] if CheckM1 skipped

    output:
    tuple val(meta), path('representatives/*.fasta'), emit: representatives, optional: true
    tuple val(meta), path('cluster_definition.tsv'),  emit: clusters
    path 'versions.yml',                              emit: versions

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
    tuple val(meta), path("*.bam"), optional: true,     emit: bams
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

/*
 * Like COVERM_GENOME, but each sample's reads are mapped to ITS OWN per-unit genome
 * set (paired via a channel join) rather than one genome set broadcast to every
 * sample. Used by the within-sample/within-group dereplication read mapping.
 */
process COVERM_GENOME_PAIRED {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads), path(genomes, stageAs: 'genomes/*')

    output:
    tuple val(meta), path("${meta.id}_abundances.tsv"), emit: abundance
    tuple val(meta), path("*.bam"), optional: true,     emit: bams
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
