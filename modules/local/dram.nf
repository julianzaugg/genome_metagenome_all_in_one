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
    # The prepared DB dir carries no CONFIG, so DRAM would otherwise fall back to
    # its packaged default (all paths unset) and annotate against nothing. Build a
    # config pointing at the staged DB files and point DRAM at it. See
    # bin/build_dram_config.py for the database selection (uniref90 excluded).
    build_dram_config.py ${dram_db} LOCAL_DRAM_CONFIG.json
    export DRAM_CONFIG_LOCATION="\$PWD/LOCAL_DRAM_CONFIG.json"

    dram.sh annotate_genes ${args} \\
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
    # See DRAM_ANNOTATE: the staged DB carries no CONFIG, so build one that
    # points at the staged DB files (uniref90 excluded) before annotating.
    build_dram_config.py ${dram_db} LOCAL_DRAM_CONFIG.json
    export DRAM_CONFIG_LOCATION="\$PWD/LOCAL_DRAM_CONFIG.json"

    dram.sh annotate_genes ${args} \\
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
    # Build the same DB config used for annotation so distill can find the
    # distillation sheets (genome_summary_form, etc_module_database, …), which
    # the staged DB has but its (absent) CONFIG never points DRAM at.
    build_dram_config.py ${dram_db} LOCAL_DRAM_CONFIG.json
    export DRAM_CONFIG_LOCATION="\$PWD/LOCAL_DRAM_CONFIG.json"

    dram.sh distill ${args} \\
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
