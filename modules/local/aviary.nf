/*
 * Aviary — metagenome bin recovery (wraps a Snakemake pipeline).
 * Mirrors run_aviary.sh `aviary recover`. Coarse-grained: the whole recover run
 * is one process (its internal steps are opaque to Nextflow -resume).
 *
 * NOTE: provide the Aviary container (complex deps). DB locations are passed as
 * flags; CONDA_ENV_PATH/env handling is expected to live inside the image.
 */

process AVIARY_RECOVER {
    tag   { meta.id }
    label 'process_maxmem'
    label 'process_long'

    input:
    tuple val(meta), path(assembly), path(reads)
    path(gtdb_db)
    path(checkm2_db)
    path(eggnog_db)

    output:
    tuple val(meta), path("${meta.id}/bins/**/*.fasta"), emit: bins
    tuple val(meta), path("${meta.id}"),                 emit: recover_dir
    path 'versions.yml',                                 emit: versions

    script:
    def args   = task.ext.args ?: '--skip-singlem'
    def reads_arg = meta.single_end ? "--longreads ${reads} --long-read-type ont" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    export GTDBTK_DATA_PATH=${gtdb_db}
    export CHECKM2_DATA_PATH=${checkm2_db}
    export EGGNOG_DATA_DIR=${eggnog_db}

    aviary recover ${args} \\
        --assembly ${assembly} \\
        ${reads_arg} \\
        --output ${meta.id} \\
        --max_threads ${task.cpus} \\
        --pplacer_threads ${Math.max(1, (task.cpus as int).intdiv(2))} \\
        --max_memory ${task.memory.toGiga()} \\
        --gtdb_path ${gtdb_db} \\
        --checkm2-db-path ${checkm2_db} \\
        --eggnog-db-path ${eggnog_db}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aviary: \$(aviary --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/bins/final_bins
    echo ">contig1" > ${meta.id}/bins/final_bins/${meta.id}_bin.1.fasta
    echo "ACGTACGT" >> ${meta.id}/bins/final_bins/${meta.id}_bin.1.fasta
    echo '"${task.process}": {aviary: stub}' > versions.yml
    """
}
