/*
 * Sylph — fast read-level taxonomic profiling.
 * Mirrors run_sylph.sh: per-sample sketch, then a combined profile across samples.
 */

process SYLPH_SKETCH {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.sylsp"), emit: sketch
    path 'versions.yml',                       emit: versions

    script:
    def args  = task.ext.args ?: ''
    def input = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph sketch ${args} -t ${task.cpus} ${input} -S ${meta.id}
    # sylph names paired sketches <sample>.paired.sylsp / single .sylsp — normalise
    mv *.sylsp ${meta.id}.sylsp 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>&1 | sed 's/sylph //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sylsp
    echo '"${task.process}": {sylph: stub}' > versions.yml
    """
}

process SYLPH_PROFILE {
    label 'process_medium'

    input:
    path(sketches, stageAs: 'sketch_*')
    path(dbs)     // one or more .syldb files (chosen via params.sylph_db glob)

    output:
    path 'sylph_profile.tsv', emit: profile
    path 'versions.yml',      emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    sylph profile ${args} -t ${task.cpus} ${dbs} ${sketches} > sylph_profile.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph: \$(sylph --version 2>&1 | sed 's/sylph //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "Sample_file\\tGenome_file\\tTaxonomic_abundance" > sylph_profile.tsv
    echo '"${task.process}": {sylph: stub}' > versions.yml
    """
}

process SYLPH_TAX {
    // Attach taxonomy to a sylph profile and merge per-sample tables.
    // Mirrors run_sylph.sh: sylph-tax taxprof -> sylph-tax merge (rel + seq abundance).
    label 'process_low'

    input:
    path(profile)
    path(metadata)   // one or more *_metadata.tsv.gz matching the .syldb files used

    output:
    path 'merged_relative_abundance.tsv', emit: relative
    path 'merged_sequence_abundance.tsv', emit: sequence
    path 'taxprof',                       emit: per_sample
    path 'versions.yml',                  emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p taxprof
    sylph-tax taxprof ${profile} -o taxprof/ -t ${metadata} ${args}
    sylph-tax merge taxprof/*.sylphmpa --column relative_abundance -o merged_relative_abundance.tsv
    sylph-tax merge taxprof/*.sylphmpa --column sequence_abundance -o merged_sequence_abundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph-tax: \$(sylph-tax --version 2>&1 | sed 's/sylph-tax //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p taxprof
    touch taxprof/sample.sylphmpa merged_relative_abundance.tsv merged_sequence_abundance.tsv
    echo '"${task.process}": {sylph-tax: stub}' > versions.yml
    """
}
