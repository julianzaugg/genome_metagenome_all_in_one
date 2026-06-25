/*
 * Gene catalogue: collect predicted proteins -> keep complete genes + namespace
 * headers by sample -> cluster with CD-HIT -> tabulate cluster membership.
 * Mirrors run_cdhit_scaffolds.sh (complete genes = pyrodigal partial=00) and the
 * downstream tabulate_cdhit.py helper. Shared by both metagenome tracks.
 */

process CATALOGUE_PREP {
    label 'process_low'

    input:
    path(faas, stageAs: 'faa/*')   // per-sample protein fastas (named <sample>.faa)

    output:
    path 'all_complete_proteins.faa', emit: proteins
    path 'versions.yml',              emit: versions

    script:
    """
    : > all_complete_proteins.faa
    for f in faa/*.faa; do
        name=\$(basename "\$f" .faa)
        # keep only complete genes (pyrodigal 'partial=00') and namespace by sample
        awk -v s="\$name" '
            /^>/ { keep = (\$0 ~ /partial=00/); if (keep) { sub(/^>/, ">" s "___"); print \$1 } next }
            { if (keep) print }
        ' "\$f" >> all_complete_proteins.faa
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/[^0-9.]*//;s/,.*//')
    END_VERSIONS
    """

    stub:
    """
    cat faa/*.faa > all_complete_proteins.faa 2>/dev/null || echo ">x___1" > all_complete_proteins.faa
    echo '"${task.process}": {awk: stub}' > versions.yml
    """
}

process CDHIT {
    label 'process_high'

    input:
    path(proteins)

    output:
    path 'gene_catalogue.faa',       emit: catalogue
    path 'gene_catalogue.faa.clstr', emit: clusters
    path 'versions.yml',             emit: versions

    script:
    def args = task.ext.args ?: '-c 1.0 -aS 0.9 -d 0 -M 0 -g 1'
    """
    cd-hit ${args} -T ${task.cpus} -i ${proteins} -o gene_catalogue.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit -h 2>&1 | grep -m1 -i version | sed 's/.*version //;s/ .*//' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    cp ${proteins} gene_catalogue.faa
    echo ">Cluster 0" > gene_catalogue.faa.clstr
    echo '"${task.process}": {cdhit: stub}' > versions.yml
    """
}

process CATALOGUE_TABULATE {
    label 'process_single'

    input:
    path(clstr)

    output:
    path 'gene_catalogue_membership.tsv', emit: membership
    path 'versions.yml',                  emit: versions

    script:
    """
    tabulate_cdhit.py -i ${clstr} > gene_catalogue_membership.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "cluster\\trepresentative\\tmember" > gene_catalogue_membership.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
