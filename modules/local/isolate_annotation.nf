/*
 * Isolate annotation and typing wrappers.
 */

process BAKTA_BAKTA {
    tag   { meta.id }
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(fasta)
    path bakta_db
    path reference_proteins

    output:
    tuple val(meta), path("${meta.id}/${meta.id}.gff3"), emit: gff
    tuple val(meta), path("${meta.id}/${meta.id}.faa"),  emit: faa
    tuple val(meta), path("${meta.id}/${meta.id}.ffn"),  emit: ffn,  optional: true
    tuple val(meta), path("${meta.id}/${meta.id}.json"), emit: json
    tuple val(meta), path("${meta.id}"),                 emit: outdir
    path 'versions.yml',                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    def proteins = reference_proteins ? "--proteins ${reference_proteins}" : ''
    """
    bakta ${args} \\
        --db ${bakta_db} \\
        --threads ${task.cpus} \\
        --prefix ${meta.id} \\
        --output ${meta.id} \\
        ${proteins} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(bakta --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    echo "##gff-version 3" > ${meta.id}/${meta.id}.gff3
    echo ">${meta.id}_cds1" > ${meta.id}/${meta.id}.faa
    echo "MAAA" >> ${meta.id}/${meta.id}.faa
    echo ">${meta.id}_cds1" > ${meta.id}/${meta.id}.ffn
    echo "ATGGCT" >> ${meta.id}/${meta.id}.ffn
    echo '{"sample":"${meta.id}"}' > ${meta.id}/${meta.id}.json
    echo '"${task.process}": {bakta: stub}' > versions.yml
    """
}

process BAKTA_STATS {
    label 'process_low'

    input:
    path jsons

    output:
    path 'bakta_annotation_stats.tsv', emit: stats
    path 'versions.yml',               emit: versions

    script:
    """
    collect-annotation-stats.py --prefix bakta_annotation_stats --output . ${jsons}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "sample\\tn_features" > bakta_annotation_stats.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

process MLST {
    label 'process_medium'

    input:
    path fastas

    output:
    path 'mlst_summary.csv', emit: summary
    path 'novel_mlst.fasta', emit: novel, optional: true
    path 'versions.yml',     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mlst ${args} --threads ${task.cpus} --novel novel_mlst.fasta --csv ${fastas} > mlst_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$(mlst --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo "FILE,SCHEME,ST" > mlst_summary.csv
    touch novel_mlst.fasta
    echo '"${task.process}": {mlst: stub}' > versions.yml
    """
}

process AMRFINDERPLUS_RUN {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(faa), path(gff)
    path amrfinder_db

    output:
    tuple val(meta), path("${meta.id}_amrfinder.tsv"), emit: results
    tuple val(meta), path("${meta.id}.amrfinder.gff"), emit: gff
    path 'versions.yml',                              emit: versions

    script:
    def args = task.ext.args ?: '--plus'
    def db = amrfinder_db ? "--database ${amrfinder_db}" : ''
    """
    grep -P "\\tCDS\\t" ${gff} | sed "s/Name=/OtherName=/g; s/ID=/Name=/g" > ${meta.id}.amrfinder.gff || true
    amrfinder ${args} \\
        --protein ${faa} \\
        --gff ${meta.id}.amrfinder.gff \\
        --threads ${task.cpus} \\
        ${db} > ${meta.id}_amrfinder.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinder: \$(amrfinder --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo -e "Protein identifier\\tGene symbol" > ${meta.id}_amrfinder.tsv
    echo "##gff-version 3" > ${meta.id}.amrfinder.gff
    echo '"${task.process}": {amrfinder: stub}' > versions.yml
    """
}

process AMRFINDERPLUS_COLLATE {
    label 'process_low'

    input:
    path tables

    output:
    path 'amrfinder_all.tsv', emit: summary
    path 'versions.yml',      emit: versions

    script:
    """
    first=1
    : > amrfinder_all.tsv
    for t in ${tables}; do
        sample=\$(basename "\$t" _amrfinder.tsv)
        if [[ "\$first" -eq 1 ]]; then
            head -n 1 "\$t" | sed 's/^/Sample\\t/' > amrfinder_all.tsv
            first=0
        fi
        tail -n +2 "\$t" | sed "s/^/\${sample}\\t/" >> amrfinder_all.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/GNU bash, version //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "Sample\\tProtein identifier\\tGene symbol" > amrfinder_all.tsv
    echo '"${task.process}": {bash: stub}' > versions.yml
    """
}

process ISESCAN {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}"), emit: outdir
    tuple val(meta), path("${meta.id}/*.tsv"), emit: tables, optional: true
    path 'versions.yml',                 emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${meta.id}
    isescan.py ${args} \\
        --seqfile ${fasta} \\
        --output ${meta.id} \\
        --nthread ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isescan: \$(isescan.py --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    echo -e "seqID\\tisBegin\\tisEnd" > ${meta.id}/${meta.id}.tsv
    echo '"${task.process}": {isescan: stub}' > versions.yml
    """
}

process ISESCAN_COLLATE {
    label 'process_low'

    input:
    path tables

    output:
    path 'ISEScan_summary.tsv', emit: summary
    path 'versions.yml',        emit: versions

    script:
    """
    first=1
    : > ISEScan_summary.tsv
    for t in ${tables}; do
        genome=\$(dirname "\$t" | xargs basename)
        if [[ "\$first" -eq 1 ]]; then
            head -n 1 "\$t" | sed 's/^/Genome_ID\\t/' > ISEScan_summary.tsv
            first=0
        fi
        tail -n +2 "\$t" | sed "s/^/\${genome}\\t/" >> ISEScan_summary.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/GNU bash, version //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "Genome_ID\\tseqID\\tisBegin\\tisEnd" > ISEScan_summary.tsv
    echo '"${task.process}": {bash: stub}' > versions.yml
    """
}
