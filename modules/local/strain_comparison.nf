/*
 * Strain-level comparison of samples against one shared reference genome set.
 *
 * Answers "are two samples recovering the SAME strain, or different strains of the
 * same species?" — which 95% ANI dereplication cannot: that only establishes a shared
 * species. Every sample's reads are mapped to ONE cross-sample reference set (the
 * dereplicated MAGs), then per-sample allele profiles are compared pairwise.
 *
 *   STRAIN_GENOME_FILTER  — MIMAG-style completeness/contamination filter on the
 *                           reference set (stricter than the pipeline HQ filter;
 *                           contamination is what drives false "different strain" calls)
 *   STRAIN_REFERENCE_PREP — combined FASTA + scaffold-to-bin (.stb), contig headers
 *                           prefixed with the bin name for cross-sample uniqueness
 *
 *   inStrain (short reads only):
 *     BOWTIE2_STRAIN_BUILD -> BOWTIE2_STRAIN_ALIGN -> INSTRAIN_PROFILE
 *       -> INSTRAIN_COMPARE -> INSTRAIN_SUMMARISE
 *
 *   TRACS (short OR long reads):
 *     TRACS_BUILD_DB -> TRACS_ALIGN -> TRACS_COMBINE -> TRACS_DISTANCE -> TRACS_CLUSTER
 */

/*
 * Restrict the reference set to genomes good enough for strain calling.
 *
 * The pipeline's HQ filter (completeness - 3*contamination >= 50) is tuned for
 * abundance estimation and admits e.g. a 95% complete / 15% contaminated bin. That is
 * too loose here: foreign contigs in a bin collect reads from unrelated populations,
 * and when those populations differ between samples they generate spurious SNPs that
 * push pairs toward a false "different strain" call. Independent completeness and
 * contamination thresholds (MIMAG high-quality by default) are the right control.
 *
 * Genomes absent from the CheckM report(s) are KEPT — this preserves the existing
 * "external reference genomes bypass the quality filter" behaviour of
 * COVERM_CLUSTER_HQ_REF when strain_genome_source = hq_ref_representatives.
 */
process STRAIN_GENOME_FILTER {
    label 'process_single'

    input:
    path(genomes, stageAs: 'genomes/*')
    path(checkm2_report)   // [] if CheckM2 skipped
    path(checkm1_report)   // [] if CheckM1 skipped

    output:
    path 'reference_genomes/*.fasta',    emit: genomes, optional: true
    path 'strain_reference_genomes.tsv', emit: report
    path 'versions.yml',                 emit: versions

    script:
    def min_completeness  = task.ext.min_completeness  ?: 90
    def max_contamination = task.ext.max_contamination ?: 5
    def hq_source         = task.ext.hq_source ?: 'both'
    """
    # Collect (genome, completeness, contamination, source) from the selected report(s).
    # Column positions differ between CheckM1's tab table and CheckM2's quality report,
    # so find them by header name -- the same idiom the CoverM clusterers use.
    : > qc_values.tsv
    read_qc() {
        [ -s "\$1" ] || return 0
        awk -F '\\t' -v src="\$2" '
            NR==1 { for (i=1; i<=NF; i++) { if (\$i=="Completeness") cc=i; if (\$i=="Contamination") ct=i;
                                            if (\$i=="Name" || \$i=="Bin Id" || \$i=="genome") id=i } next }
            (cc && ct && id) { print \$id "\\t" \$cc "\\t" \$ct "\\t" src }
        ' "\$1" >> qc_values.tsv
    }
    case "${hq_source}" in
        both)    read_qc "${checkm2_report}" checkm2; read_qc "${checkm1_report}" checkm1 ;;
        checkm2) read_qc "${checkm2_report}" checkm2 ;;
        checkm1) read_qc "${checkm1_report}" checkm1 ;;
    esac

    # Genome names, from the staged fasta basenames.
    : > genome_names.txt
    for f in genomes/*.fasta genomes/*.fa genomes/*.fna; do
        [ -e "\$f" ] || continue
        b=\$(basename "\$f"); printf '%s\\n' "\${b%.*}" >> genome_names.txt
    done
    sort -u genome_names.txt -o genome_names.txt

    # One verdict per genome, in a single pass. A genome passes if ANY selected report
    # passes it (matching hq_quality_source semantics elsewhere in the pipeline). A
    # genome present in no report is KEPT and flagged 'unscored', so user-supplied
    # reference genomes -- which never appear in the MAG CheckM reports -- are not
    # silently dropped.
    awk -F '\\t' -v mc=${min_completeness} -v xc=${max_contamination} '
        FNR==NR {
            g=\$1; ok = (\$2 >= mc && \$3 <= xc)
            if (!(g in seen) || (ok && !pass[g])) { comp[g]=\$2; cont[g]=\$3; src[g]=\$4; pass[g]=ok }
            seen[g]=1
            next
        }
        {
            g=\$1
            if (!(g in seen))    { print g "\\tNA\\tNA\\tnone\\tunscored\\tnot_in_checkm_report" }
            else if (pass[g])    { print g "\\t" comp[g] "\\t" cont[g] "\\t" src[g] "\\tkept\\tpassed" }
            else                 { print g "\\t" comp[g] "\\t" cont[g] "\\t" src[g] "\\tdropped\\tbelow_thresholds" }
        }
    ' qc_values.tsv genome_names.txt > verdicts.tsv

    printf 'genome\\tcompleteness\\tcontamination\\tsource\\tstatus\\treason\\n' > strain_reference_genomes.tsv
    cat verdicts.tsv >> strain_reference_genomes.tsv

    mkdir -p reference_genomes
    awk -F '\\t' '\$5=="kept" || \$5=="unscored" { print \$1 }' verdicts.tsv | while read -r genome; do
        for f in "genomes/\${genome}.fasta" "genomes/\${genome}.fa" "genomes/\${genome}.fna"; do
            [ -e "\$f" ] && cp "\$f" "reference_genomes/\${genome}.fasta" && break
        done
    done

    n_kept=\$(awk -F '\\t' '\$5=="kept"' verdicts.tsv | wc -l)
    n_unscored=\$(awk -F '\\t' '\$5=="unscored"' verdicts.tsv | wc -l)
    n_dropped=\$(awk -F '\\t' '\$5=="dropped"' verdicts.tsv | wc -l)
    echo "strain reference set: \${n_kept} passed (completeness >= ${min_completeness}, contamination <= ${max_contamination}; source ${hq_source}), \${n_dropped} dropped, \${n_unscored} kept unscored" >&2
    if [ "\$((n_kept + n_unscored))" -eq 0 ]; then
        echo "WARNING: no genomes passed the strain reference quality filter -- strain comparison will be skipped. See strain_reference_genomes.tsv; consider relaxing --strain_min_completeness / --strain_max_contamination." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p reference_genomes
    cp genomes/* reference_genomes/ 2>/dev/null || (printf '>c\\nACGT\\n' > reference_genomes/rep.1.fasta)
    printf 'genome\\tcompleteness\\tcontamination\\tsource\\tstatus\\treason\\n' > strain_reference_genomes.tsv
    printf 'rep.1\\t99.0\\t0.5\\tcheckm2\\tkept\\tpassed\\n' >> strain_reference_genomes.tsv
    echo '"${task.process}": {awk: stub}' > versions.yml
    """
}

/*
 * Combined reference FASTA + .stb for inStrain.
 *
 * Contig headers are prefixed with the bin name. This is NOT cosmetic: bins come from
 * per-sample metaSPAdes assemblies whose contigs are named NODE_<n>_length_..._cov_...,
 * so identical names can recur across independently assembled samples. inStrain keys
 * everything on scaffold name, so a collision would silently merge two unrelated
 * contigs. (This is also why inStrain's bundled parse_stb.py is not used -- it builds
 * the .stb but does not rename contigs.)
 */
process STRAIN_REFERENCE_PREP {
    label 'process_single'

    input:
    path(genomes, stageAs: 'genomes/*')

    output:
    path 'strain_reference.fasta', emit: fasta
    path 'strain_reference.stb',   emit: stb
    path 'versions.yml',           emit: versions

    script:
    """
    : > strain_reference.fasta
    : > strain_reference.stb

    for f in genomes/*.fasta; do
        [ -e "\$f" ] || continue
        genome=\$(basename "\$f" .fasta)

        awk -v b="\$genome" '/^>/ { print ">" b "__" substr(\$0,2); next } { print }' \\
            "\$f" >> strain_reference.fasta

        awk -v b="\$genome" '/^>/ { h=substr(\$0,2); split(h,a," "); print b "__" a[1] "\\t" b }' \\
            "\$f" >> strain_reference.stb
    done

    echo "combined reference: \$(grep -c '^>' strain_reference.fasta) scaffolds across \$(cut -f2 strain_reference.stb | sort -u | wc -l) genomes" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf '>rep.1__c1\\nACGT\\n' > strain_reference.fasta
    printf 'rep.1__c1\\trep.1\\n'  > strain_reference.stb
    echo '"${task.process}": {awk: stub}' > versions.yml
    """
}

process BOWTIE2_STRAIN_BUILD {
    label 'process_high'

    input:
    path(reference)

    output:
    path 'bt2_index',    emit: index
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p bt2_index
    bowtie2-build ${args} --threads ${task.cpus} ${reference} bt2_index/strain_reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bt2_index
    touch bt2_index/strain_reference.1.bt2
    echo '"${task.process}": {bowtie2: stub}' > versions.yml
    """
}

/*
 * Emits a COORDINATE-SORTED, INDEXED bam. The .bai matters beyond convenience:
 * inStrain's prepare_bam_fie() only skips its own `samtools sort/index` shell-outs when
 * an index is already present, and the inStrain biocontainer ships pysam but NOT the
 * samtools binary. Dropping the .bai would make INSTRAIN_PROFILE fail at run time.
 */
process BOWTIE2_STRAIN_ALIGN {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path 'versions.yml',                                                 emit: versions

    script:
    def args      = task.ext.args ?: ''
    def reads_arg = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bowtie2 ${args} \\
        -p ${task.cpus} \\
        -x ${index}/strain_reference \\
        ${reads_arg} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index -@ ${task.cpus} ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/^.*version //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    echo '"${task.process}": {bowtie2: stub, samtools: stub}' > versions.yml
    """
}

process INSTRAIN_PROFILE {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(stb)

    output:
    tuple val(meta), path("${meta.id}.IS"), emit: profile
    path 'versions.yml',                    emit: versions

    script:
    def args = task.ext.args ?: '--database_mode --skip_plot_generation'
    """
    inStrain profile ${bam} ${reference} ${args} \\
        -o ${meta.id}.IS \\
        -p ${task.cpus} \\
        -s ${stb}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(inStrain --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}.IS/output
    printf 'genome\\tcoverage\\tbreadth\\n' > ${meta.id}.IS/output/${meta.id}.IS_genome_info.tsv
    echo '"${task.process}": {instrain: stub}' > versions.yml
    """
}

process INSTRAIN_COMPARE {
    label 'process_high'

    input:
    path(profiles, stageAs: 'profiles/*')
    path(stb)

    output:
    path 'strain_compare.IS',      emit: compare
    path 'excluded_profiles.tsv',  emit: excluded
    path 'versions.yml',           emit: versions

    script:
    def args = task.ext.args ?: '--database_mode --skip_plot_generation -ani 0.99999 -cov 0.5'
    """
    # --database_mode requires genome-level information in EVERY profile, and
    # aborts the entire comparison if any one lacks it (compare_utils.py:
    # find_relevant_scaffolds raises when SNVprofile.get('genome_level_info')
    # returns None). A sample with almost nothing mapping to the reference set
    # never gets that information -- database mode pins --min_genome_coverage 1,
    # so if no genome reaches 1x there is nothing to profile at genome level.
    #
    # Losing every other sample's comparison to the weakest sample in the set is
    # not a useful failure mode, so drop such profiles explicitly and record
    # which ones. SNVprofile.get() looks the name up in raw_data/attributes.tsv,
    # so its presence there is exactly the condition inStrain itself tests.
    printf 'sample\\treason\\n' > excluded_profiles.tsv
    usable=""
    n_usable=0
    for prof in profiles/*; do
        [ -d "\$prof" ] || continue
        name=\$(basename "\$prof" .IS)
        if [ -s "\$prof/raw_data/attributes.tsv" ] && \\
           awk -F '\\t' '\$1=="genome_level_info" { found=1 } END { exit !found }' "\$prof/raw_data/attributes.tsv"; then
            usable="\$usable \$prof"
            n_usable=\$((n_usable + 1))
        else
            printf '%s\\tno_genome_level_info\\n' "\$name" >> excluded_profiles.tsv
            echo "WARNING: excluding '\$name' from inStrain compare -- the profile has no genome-level information, i.e. no genome reached the 1x coverage floor that --database_mode applies. Usually means very little of this sample's reads map to the reference genome set; check its mapping rate in the read-stat report." >&2
        fi
    done

    if [ "\$n_usable" -lt 2 ]; then
        echo "ERROR: only \$n_usable profile(s) carry genome-level information, but a pairwise comparison needs at least 2. See excluded_profiles.tsv. Either too little of these samples' reads map to the reference genome set, or the set is too small/too strictly filtered (--strain_min_completeness / --strain_max_contamination)." >&2
        exit 1
    fi

    echo "inStrain compare: using \$n_usable profile(s); \$(( \$(wc -l < excluded_profiles.tsv) - 1 )) excluded" >&2

    inStrain compare -i \$usable ${args} \\
        -o strain_compare.IS \\
        -p ${task.cpus} \\
        -s ${stb}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instrain: \$(inStrain --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p strain_compare.IS/output
    printf 'genome\\tname1\\tname2\\tpopANI\\tconANI\\tpercent_genome_compared\\tcoverage_overlap\\n' \\
        > strain_compare.IS/output/strain_compare.IS_genomeWide_compare.tsv
    printf 'genome\\tcluster\\tsample\\n' \\
        > strain_compare.IS/output/strain_compare.IS_strain_clusters.tsv
    printf 'sample\\treason\\n' > excluded_profiles.tsv
    echo '"${task.process}": {instrain: stub}' > versions.yml
    """
}

process INSTRAIN_SUMMARISE {
    label 'process_single'

    input:
    path(compare_dir)

    output:
    path 'summary/*',    emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    instrain_summarise.py --compare-dir ${compare_dir} --outdir summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p summary
    printf 'genome\\tsample_a\\tsample_b\\tpopANI\\tsame_strain\\n' > summary/strain_sharing_summary.tsv
    printf 'sample_a\\tsample_b\\tn_genomes_compared\\tn_same_strain\\n' > summary/strain_sharing_counts.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

/*
 * TRACS reference database built from OUR OWN genomes -- the documented custom-database
 * path, not the GTDB/sourmash one. build-db writes a zip embedding both the sourmash
 * index and a gzipped copy of every genome, so `tracs align` extracts references
 * straight out of it: no --refseqs, no Genbank download, no network access at run time.
 *
 * Note build-db appends ".zip" to -o, so the argument must NOT already end in .zip.
 */
process TRACS_BUILD_DB {
    label 'process_high'

    input:
    path(genomes, stageAs: 'genomes/*')

    output:
    path 'strain_db.zip', emit: db
    path 'versions.yml',  emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    tracs build-db ${args} -i genomes/*.fasta -o strain_db -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tracs: \$(tracs --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch strain_db.zip
    echo '"${task.process}": {tracs: stub}' > versions.yml
    """
}

process TRACS_ALIGN {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_align"), emit: alignment
    path 'versions.yml',                       emit: versions

    script:
    def args      = task.ext.args ?: '--minimap_preset sr --keep-all'
    def reads_arg = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    mkdir -p ${meta.id}_align
    tracs align ${args} \\
        -i ${reads_arg} \\
        --database ${db} \\
        -o ${meta.id}_align \\
        -p ${meta.id} \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tracs: \$(tracs --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_align
    printf '>${meta.id}\\nACGT\\n' > ${meta.id}_align/posterior_counts_ref_rep.1.fasta
    echo '"${task.process}": {tracs: stub}' > versions.yml
    """
}

/*
 * `tracs align` writes one posterior_counts_ref_<ref>.fasta per sample; `combine`
 * transposes that into one MSA per reference genome across all samples, which is what
 * `distance` consumes. distance cannot read raw align output.
 */
process TRACS_COMBINE {
    label 'process_medium'

    input:
    path(alignments, stageAs: 'alignments/*')

    output:
    path 'combined',     emit: combined
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    tracs combine ${args} -i alignments/* -o combined -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tracs: \$(tracs --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p combined
    printf '>sample1\\nACGT\\n>sample2\\nACGT\\n' > combined/rep.1.fasta
    echo '"${task.process}": {tracs: stub}' > versions.yml
    """
}

process TRACS_DISTANCE {
    label 'process_high'

    input:
    path(combined)

    output:
    path 'transmission_distances.csv', emit: distances
    path 'versions.yml',               emit: versions

    script:
    def args = task.ext.args ?: '--filter -D 1000'
    """
    tracs distance ${args} \\
        --msa ${combined}/*.fasta \\
        -o transmission_distances.csv \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tracs: \$(tracs --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf 'sampleA,sampleB,SNP distance,filtered SNP distance,sites considered\\n' > transmission_distances.csv
    echo '"${task.process}": {tracs: stub}' > versions.yml
    """
}

process TRACS_CLUSTER {
    label 'process_single'

    input:
    path(distances)

    output:
    path 'strain_clusters.csv', emit: clusters
    path 'versions.yml',        emit: versions

    script:
    def args = task.ext.args ?: '-D filter -c 10'
    """
    tracs cluster ${args} -d ${distances} -o strain_clusters.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tracs: \$(tracs --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf 'sample,cluster\\n' > strain_clusters.csv
    echo '"${task.process}": {tracs: stub}' > versions.yml
    """
}
