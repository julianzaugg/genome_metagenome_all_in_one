/*
 * Marker-gene phylogenetic tree of MAGs + selected GTDB references.
 *   MARKER_TREE_PREP        — select references + assemble per-domain MSA
 *                             (bin/build_marker_tree_msa.py; ports the example
 *                             get_{closest,related}_reference_genomes.sh logic)
 *   MARKER_TREE_VERYFASTTREE — fast approximate-ML tree (default builder)
 *   MARKER_TREE_IQTREE       — maximum-likelihood tree (alternative builder)
 */

process MARKER_TREE_PREP {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(gtdbtk_dir)
    path(checkm2_report)
    path(genomes, stageAs: 'genomes/*')
    path(accessions)
    val(opts)

    output:
    path '*.marker_msa.fasta',        emit: msa,        optional: true
    path '*.reference_genomes.tsv',   emit: references, optional: true
    path '*.closest_references.tsv',  emit: closest,    optional: true
    path 'versions.yml',              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def checkm2_arg = checkm2_report ? "--checkm2 ${checkm2_report}" : ''
    def acc_arg     = accessions     ? "--accessions ${accessions}"  : ''
    """
    build_marker_tree_msa.py \\
        --gtdbtk-dir ${gtdbtk_dir} \\
        --genomes-dir genomes \\
        ${checkm2_arg} \\
        ${acc_arg} \\
        ${opts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        ete3: \$(python -c 'import ete3; print(ete3.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    cat <<-MSA > bac120.marker_msa.fasta
    >bin.1 d__Bacteria
    MKKLLAAALLAA
    >bin.2 d__Bacteria
    MKKLLAAALLAB
    >GB_GCA_000000001.1 d__Bacteria
    MKKLLAAALLAC
    MSA
    echo -e "bin.1\\tGB_GCA_000000001.1" > bac120.closest_references.tsv
    echo -e "GB_GCA_000000001.1\\td__Bacteria" > bac120.reference_genomes.tsv
    echo '"${task.process}": {python: stub, ete3: stub}' > versions.yml
    """
}

process MARKER_TREE_VERYFASTTREE {
    tag   { domain }
    label 'process_high'

    input:
    tuple val(domain), path(msa)

    output:
    tuple val(domain), path("${domain}.treefile"), emit: tree
    path 'versions.yml',                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vft=\$(command -v VeryFastTree || command -v veryfasttree)
    "\$vft" -threads ${task.cpus} ${args} ${msa} > ${domain}.treefile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        veryfasttree: \$("\$vft" --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo "(bin.1,bin.2,GB_GCA_000000001.1);" > ${domain}.treefile
    echo '"${task.process}": {veryfasttree: stub}' > versions.yml
    """
}

process MARKER_TREE_IQTREE {
    tag   { domain }
    label 'process_high'

    input:
    tuple val(domain), path(msa)

    output:
    tuple val(domain), path("${domain}.treefile"), emit: tree, optional: true
    tuple val(domain), path("${domain}.*"),        emit: outputs
    path 'versions.yml',                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    iqtree -s ${msa} -pre ${domain} -nt ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo "(bin.1,bin.2,GB_GCA_000000001.1);" > ${domain}.treefile
    echo '"${task.process}": {iqtree: stub}' > versions.yml
    """
}
