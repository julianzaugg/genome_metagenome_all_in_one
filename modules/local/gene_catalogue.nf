/*
 * Gene catalogue: collect predicted proteins + nucleotides -> keep complete genes
 * (pyrodigal partial=00) namespaced by sample -> cluster with CD-HIT at 100% and
 * 90% -> extract nucleotide CDS for the 100% reps -> tabulate cluster membership.
 * Mirrors run_cdhit_scaffolds.sh (mfqe replaced by seqkit). Shared by both
 * metagenome tracks.
 */

process CATALOGUE_PREP {
    label 'process_low'

    input:
    path(faas, stageAs: 'faa/*')   // per-sample protein fastas (<sample>.faa)
    path(fnas, stageAs: 'fna/*')   // per-sample nucleotide fastas (<sample>.fna)

    output:
    path 'all_complete_proteins.faa',    emit: proteins
    path 'all_complete_nucleotides.fna', emit: nucleotides
    path 'versions.yml',                 emit: versions

    script:
    """
    : > all_complete_proteins.faa
    : > all_complete_nucleotides.fna
    for f in faa/*.faa; do
        name=\$(basename "\$f" .faa)
        # keep complete genes (header contains partial=00) and namespace by sample
        seqkit grep -n -r -p 'partial=00' "\$f"            | sed "s/^>/>\${name}___/" >> all_complete_proteins.faa
        seqkit grep -n -r -p 'partial=00' "fna/\${name}.fna" | sed "s/^>/>\${name}___/" >> all_complete_nucleotides.fna
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    cat faa/*.faa > all_complete_proteins.faa 2>/dev/null || echo ">x___1" > all_complete_proteins.faa
    cat fna/*.fna > all_complete_nucleotides.fna 2>/dev/null || echo ">x___1" > all_complete_nucleotides.fna
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}

process CDHIT {
    label 'process_high'

    input:
    path(proteins)
    val(identities)   // comma-separated, e.g. "1.0,0.9" — first is the primary catalogue

    output:
    path 'gene_catalogue.faa',       emit: catalogue       // primary (first identity)
    path 'gene_catalogue.faa.clstr', emit: clusters        // primary .clstr
    path 'gene_catalogue_id*.faa',   emit: catalogues      // one per identity
    path 'versions.yml',             emit: versions

    script:
    def args = task.ext.args ?: '-d 0 -M 80000 -g 1'
    """
    first=1
    for c in \$(echo "${identities}" | tr ',' ' '); do
        # cd-hit word size depends on the identity threshold
        n=\$(awk -v c="\$c" 'BEGIN{ if(c>=0.7)print 5; else if(c>=0.6)print 4; else if(c>=0.5)print 3; else print 2 }')
        pct=\$(awk -v c="\$c" 'BEGIN{ printf "%d", c*100 }')
        cd-hit ${args} -c "\$c" -n "\$n" -T ${task.cpus} -i ${proteins} -o gene_catalogue_id\${pct}.faa
        if [ "\$first" = 1 ]; then
            cp gene_catalogue_id\${pct}.faa       gene_catalogue.faa
            cp gene_catalogue_id\${pct}.faa.clstr gene_catalogue.faa.clstr
            first=0
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cdhit: \$(cd-hit -h 2>&1 | grep -m1 -i version | sed 's/.*version //;s/ .*//' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    for c in \$(echo "${identities}" | tr ',' ' '); do
        pct=\$(awk -v c="\$c" 'BEGIN{ printf "%d", c*100 }')
        cp ${proteins} gene_catalogue_id\${pct}.faa
    done
    cp ${proteins} gene_catalogue.faa
    echo ">Cluster 0" > gene_catalogue.faa.clstr
    echo '"${task.process}": {cdhit: stub}' > versions.yml
    """
}

process CATALOGUE_CDS {
    // nucleotide sequences for the 100% representative proteins (seqkit grep)
    label 'process_low'

    input:
    path(catalogue)     // gene_catalogue.faa (100% reps)
    path(nucleotides)   // all_complete_nucleotides.fna

    output:
    path 'gene_catalogue_CDS.fna', emit: cds
    path 'versions.yml',           emit: versions

    script:
    """
    grep '^>' ${catalogue} | sed 's/^>//; s/ .*//' > rep_ids.txt
    seqkit grep -f rep_ids.txt ${nucleotides} > gene_catalogue_CDS.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    touch gene_catalogue_CDS.fna
    echo '"${task.process}": {seqkit: stub}' > versions.yml
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
