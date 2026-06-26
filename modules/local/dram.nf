/*
 * DRAM — functional annotation of the gene catalogue proteins.
 * Mirrors run_dram_cdhit.sh (DRAM.py annotate_genes on a protein fasta).
 */

process DRAM_ANNOTATE {
    label 'process_high'
    label 'process_long'

    input:
    path(proteins)
    path(dram_db)

    output:
    path 'dram_annotations/annotations.tsv', emit: annotations
    path 'dram_annotations',                 emit: outdir
    path 'versions.yml',                     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # DRAM reads DB locations from its config. Some installs keep CONFIG
    # outside the database directory, so only override it for prepared DB dirs
    # that actually carry their own CONFIG file.
    if [[ -f "${dram_db}/CONFIG" ]]; then
        export DRAM_CONFIG_LOCATION="${dram_db}/CONFIG"
    fi

    DRAM.py annotate_genes ${args} \\
        --input_faa ${proteins} \\
        --output_dir dram_annotations \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: \$(DRAM.py --version 2>&1 | sed 's/DRAM version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dram_annotations
    echo -e "\\tgene_id\\tannotation" > dram_annotations/annotations.tsv
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}
