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
        dram: \$(pip show dram 2>/dev/null | grep '^Version:' | sed 's/Version: //' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dram_annotations
    echo -e "\\tgene_id\\tannotation" > dram_annotations/annotations.tsv
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}

/*
 * DRAM — per-bin functional annotation. Mirrors run_dram_bins_parallel.sh
 * (DRAM.py annotate_genes on each bin's .faa). Nextflow's natural per-item
 * parallelism replaces GNU Parallel.
 */

process DRAM_ANNOTATE_BINS {
    tag { meta.id }
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(proteins)
    path(dram_db)

    output:
    tuple val(meta), path('dram_annotations/annotations.tsv'), emit: annotations
    tuple val(meta), path('dram_annotations'),                  emit: outdir
    path 'versions.yml',                                        emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    if [[ -f "${dram_db}/CONFIG" ]]; then
        export DRAM_CONFIG_LOCATION="${dram_db}/CONFIG"
    fi

    DRAM.py annotate_genes ${args} \\
        --input_faa ${proteins} \\
        --output_dir dram_annotations \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: \$(pip show dram 2>/dev/null | grep '^Version:' | sed 's/Version: //' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dram_annotations
    echo -e "\\tgene_id\\tannotation" > dram_annotations/annotations.tsv
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}

process DRAM_DISTILL {
    tag { meta.id }
    label 'process_low'

    input:
    tuple val(meta), path(annotations)
    path(dram_db)

    output:
    tuple val(meta), path('distilled'), emit: distilled
    path 'versions.yml',                emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    if [[ -f "${dram_db}/CONFIG" ]]; then
        export DRAM_CONFIG_LOCATION="${dram_db}/CONFIG"
    fi

    DRAM.py distill ${args} \\
        --input_file ${annotations} \\
        --output_dir distilled

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: \$(pip show dram 2>/dev/null | grep '^Version:' | sed 's/Version: //' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p distilled
    touch distilled/metabolism_summary.xlsx
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}
