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
    tuple val(meta), path("${meta.id}/renamed_bins/*.fasta"),             emit: bins
    tuple val(meta), path("${meta.id}/renamed_bins/bin_contig_list.tsv"), emit: bin_contig_list
    tuple val(meta), path("${meta.id}/bins/**/*.fasta", optional: true),  emit: raw_bins
    tuple val(meta), path("${meta.id}"),                                  emit: recover_dir
    path 'versions.yml',                                                  emit: versions

    script:
    def args   = task.ext.args ?: '--skip-singlem'
    def long_read_type = task.ext.long_read_type ?: 'ont'
    def reads_arg = meta.single_end ? "--longreads ${reads} --long-read-type ${long_read_type}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def aviary_container_hint = params.aviary_container ?: "${params.container_base}/aviary_0.13.0.sif"
    """
    if ! command -v pixi >/dev/null 2>&1; then
        cat >&2 <<'EOF'
ERROR: Aviary 0.13.0 requires pixi inside the AVIARY_RECOVER runtime image.
The quay.io/biocontainers/aviary:0.13.0--pyhdfd78af_0 image contains the Aviary
CLI but not pixi, so Aviary's internal Snakemake rules fail with:
  /usr/bin/bash: line 1: pixi: command not found

Build/provide an upstream-style Aviary SIF with pixi and prebuilt Aviary
environments, then set --aviary_container if it is not at:
  ${aviary_container_hint}
EOF
        exit 127
    fi

    # Nextflow runs via 'apptainer exec' which bypasses ENTRYPOINT.
    # aviary lives in the pixi default env; add it to PATH explicitly.
    export PATH="/aviary/.pixi/envs/default/bin:\${PATH}"
    export CONDA_PREFIX="/aviary/.pixi/envs/default"

    # Resolve absolute paths — aviary's internal Snakemake changes working
    # directory, so staged relative paths (e.g. "CheckM2_database") break.
    gtdb_abs=\$(realpath ${gtdb_db})
    checkm2_abs=\$(realpath ${checkm2_db})
    eggnog_abs=\$(realpath ${eggnog_db})

    # CHECKM2DB is what checkm2 actually reads; CHECKM2_DATA_PATH is not used.
    export CHECKM2DB="\${checkm2_abs}"
    export GTDBTK_DATA_PATH="\${gtdb_abs}"
    export EGGNOG_DATA_DIR="\${eggnog_abs}"
    export SINGLEM_METAPACKAGE_PATH=\${SINGLEM_METAPACKAGE_PATH:-not_required_with_skip_singlem}
    export METABULI_DB_PATH=\${METABULI_DB_PATH:-not_required_by_gmaio_aviary_recover}
    # Redirect pixi's repodata cache off NFS (Bunya home dirs) onto local /tmp
    export PIXI_CACHE_DIR=\${PIXI_CACHE_DIR:-/tmp/pixi-cache-\${USER:-runner}}

    aviary recover ${args} \\
        --assembly ${assembly} \\
        ${reads_arg} \\
        --output ${meta.id} \\
        --max_threads ${task.cpus} \\
        --pplacer_threads ${Math.max(1, (task.cpus as int).intdiv(2))} \\
        --max_memory ${task.memory.toGiga()} \\
        --gtdb_path "\${gtdb_abs}" \\
        --checkm2-db-path "\${checkm2_abs}" \\
        --eggnog-db-path "\${eggnog_abs}"

    mkdir -p ${meta.id}/renamed_bins
    shopt -s nullglob
    for bin_file in ${meta.id}/bins/final_bins/*.{fna,fa,fasta}; do
        bin_file_basename=\$(basename "\$bin_file")
        bin_file_basename="\${bin_file_basename%.fasta}"
        bin_file_basename="\${bin_file_basename%.fna}"
        bin_file_basename="\${bin_file_basename%.fa}"
        bin_file_basename="\${bin_file_basename%.tsv}"

        binner=\$(echo "\$bin_file_basename" | awk -F "." '{print \$1}' | sed 's/_bins\$//')
        bin_number=\$(echo "\$bin_file_basename" | awk -F "." '{print \$2}')

        if [[ "\$bin_file_basename" == *"rosella_dastool_refined"* ]]; then
            bin_number="\${bin_file_basename#*refined_}"
            binner="rosella_dastool_refined"
        fi
        if [[ "\$bin_file_basename" == *"single_contig_dastool_refined"* ]]; then
            bin_number="\${bin_file_basename#*refined_}"
            binner="single_contig_dastool_refined"
        fi

        if [[ -n "\$bin_number" ]]; then
            bin_out_name="\${binner}.\${bin_number}"
        else
            bin_out_name="\${binner}"
        fi

        cp -L "\$bin_file" "${meta.id}/renamed_bins/${meta.id}.\${bin_out_name}.fasta"
    done

    : > ${meta.id}/renamed_bins/bin_contig_list.tsv
    for i in ${meta.id}/renamed_bins/*.fasta; do
        name=\$(basename "\$i" .fasta)
        grep '^>' "\$i" | sed "s/^>/\${name}\\t/" >> ${meta.id}/renamed_bins/bin_contig_list.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aviary: \$(aviary --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/bins/final_bins
    echo ">contig1" > ${meta.id}/bins/final_bins/metabat2_refined.001.fasta
    echo "ACGTACGT" >> ${meta.id}/bins/final_bins/metabat2_refined.001.fasta
    mkdir -p ${meta.id}/renamed_bins
    cp ${meta.id}/bins/final_bins/metabat2_refined.001.fasta ${meta.id}/renamed_bins/${meta.id}.metabat2_refined.001.fasta
    echo -e "${meta.id}.metabat2_refined.001\\tcontig1" > ${meta.id}/renamed_bins/bin_contig_list.tsv
    echo '"${task.process}": {aviary: stub}' > versions.yml
    """
}

process AVIARY_COLLECT_BINS {
    label 'process_low'

    input:
    path(bins, stageAs: 'input_bins/*')

    output:
    path 'all_aviary_bins/*.fasta',          emit: bins
    path 'all_aviary_bins/bin_contig_list.tsv', emit: bin_contig_list
    path 'versions.yml',                     emit: versions

    script:
    """
    mkdir -p all_aviary_bins
    cp -L input_bins/*.fasta all_aviary_bins/

    : > all_aviary_bins/bin_contig_list.tsv
    for i in all_aviary_bins/*.fasta; do
        name=\$(basename "\$i" .fasta)
        grep '^>' "\$i" | sed "s/^>/\${name}\\t/" >> all_aviary_bins/bin_contig_list.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/GNU bash, version //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p all_aviary_bins
    cp -L input_bins/*.fasta all_aviary_bins/
    : > all_aviary_bins/bin_contig_list.tsv
    for i in all_aviary_bins/*.fasta; do
        name=\$(basename "\$i" .fasta)
        grep '^>' "\$i" | sed "s/^>/\${name}\\t/" >> all_aviary_bins/bin_contig_list.tsv
    done
    echo '"${task.process}": {bash: stub}' > versions.yml
    """
}
