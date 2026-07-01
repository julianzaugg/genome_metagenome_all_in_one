/*
 * DRAM — functional annotation + distillation.
 *
 * The prepared DRAM database directory carries no CONFIG (it lived in the setup
 * host's mag_annotator package), so every process rebuilds one from the staged DB
 * with build_dram_config.py and runs DRAM through the dram.sh wrapper (which also
 * patches known container bugs). See bin/build_dram_config.py and bin/dram.sh.
 *
 * Gene catalogue: the catalogue .faa is split (DRAM_ANNOTATE_SPLIT) and annotated
 * in parallel (DRAM_ANNOTATE per chunk), then reassembled (DRAM_ANNOTATE_MERGE) —
 * ORF annotation is independent per gene, so this is exact for the hmm/pfam
 * databases (minor peptidase reciprocal-best-hit differences at chunk edges).
 * Bins: each bin is annotated separately (DRAM_ANNOTATE_BINS), then all per-bin
 * annotations are concatenated (DRAM_COMBINE_ANNOTATIONS) and distilled once
 * (DRAM_DISTILL) into a single cross-MAG summary.
 */

process DRAM_ANNOTATE_SPLIT {
    label 'process_low'

    input:
    path(proteins)

    output:
    path 'chunk.*.faa',  emit: chunks
    path 'versions.yml', emit: versions

    script:
    // Number of chunks to scatter the catalogue across (set via ext.args).
    def n = task.ext.args ?: '16'
    """
    split_fasta.py ${proteins} ${n} chunk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        split_fasta: 1
    END_VERSIONS
    """

    stub:
    """
    printf '>g1\\nMAS\\n' > chunk.001.faa
    printf '>g2\\nMAS\\n' > chunk.002.faa
    echo '"${task.process}": {split_fasta: stub}' > versions.yml
    """
}

process DRAM_ANNOTATE {
    tag { proteins.baseName }
    label 'process_high'
    label 'process_long'

    input:
    path(proteins)
    path(dram_db)

    output:
    path 'dram_annotations/annotations.tsv',      emit: annotations
    path 'dram_annotations/genes.annotated.faa',  emit: annotated_faa
    path 'versions.yml',                          emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # The prepared DB dir carries no CONFIG, so DRAM would otherwise fall back to
    # its packaged default (all paths unset) and annotate against nothing. Build a
    # config pointing at the staged DB files and point DRAM at it. See
    # bin/build_dram_config.py for the database selection (uniref90 + pfam excluded).
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
    echo -e "\\tfasta\\tko_id" > dram_annotations/annotations.tsv
    echo -e ">${proteins.baseName}_gene1\\nMASE" > dram_annotations/genes.annotated.faa
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}

process DRAM_ANNOTATE_MERGE {
    label 'process_low'

    input:
    path(annotations, stageAs: 'anno_*.tsv')
    path(faas,        stageAs: 'faa_*.faa')

    output:
    path 'dram_annotations/annotations.tsv',     emit: annotations
    path 'dram_annotations',                     emit: outdir
    path 'versions.yml',                         emit: versions

    script:
    """
    mkdir -p dram_annotations
    combine_dram_annotations.py dram_annotations/annotations.tsv --glob 'anno_*.tsv'
    cat faa_*.faa > dram_annotations/genes.annotated.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dram_annotations
    echo -e "\\tfasta\\tko_id" > dram_annotations/annotations.tsv
    touch dram_annotations/genes.annotated.faa
    echo '"${task.process}": {pandas: stub}' > versions.yml
    """
}

/*
 * DRAM — per-bin functional annotation. One task per bin; Nextflow's per-item
 * parallelism replaces the original GNU Parallel loop.
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
    # points at the staged DB files (uniref90 + pfam excluded) before annotating.
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
    echo -e "\\tfasta\\tko_id" > dram_annotations/annotations.tsv
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}

/*
 * Concatenate per-bin annotations into one table (union of columns) so a single
 * distill produces a cross-MAG summary. Each per-bin annotations.tsv already
 * carries a distinct 'fasta' column (the bin id), which distill groups on.
 */

process DRAM_COMBINE_ANNOTATIONS {
    label 'process_low'

    input:
    path(annotations, stageAs: 'anno_*.tsv')

    output:
    path 'combined_annotations.tsv', emit: annotations
    path 'versions.yml',             emit: versions

    script:
    """
    combine_dram_annotations.py combined_annotations.tsv --glob 'anno_*.tsv'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo -e "\\tfasta\\tko_id" > combined_annotations.tsv
    echo '"${task.process}": {pandas: stub}' > versions.yml
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
