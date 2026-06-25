/*
 * Reference-database downloaders. Each process uses the tool's own downloader and
 * writes into a named subdirectory of params.db_outdir (handled via publishDir).
 * Selected by the DOWNLOAD_DBS workflow from --download_db. Confirm exact
 * downloader flags against your tool versions — see docs/databases.md.
 */

process DL_GTDBTK {
    label 'process_medium'
    output:
    path 'gtdbtk', emit: db
    script:
    """
    mkdir -p gtdbtk
    download-db.sh gtdbtk 2>/dev/null || gtdbtk download-db --path gtdbtk
    """
    stub:
    "mkdir -p gtdbtk/metadata && echo VERSION_DATA=stub > gtdbtk/metadata/metadata.txt"
}

process DL_CHECKM2 {
    label 'process_low'
    output:
    path 'checkm2', emit: db
    script:
    "checkm2 database --download --path checkm2"
    stub:
    "mkdir -p checkm2"
}

process DL_GENOMAD {
    label 'process_low'
    output:
    path 'genomad', emit: db
    script:
    "genomad download-database genomad"
    stub:
    "mkdir -p genomad"
}

process DL_CHECKV {
    label 'process_low'
    output:
    path 'checkv', emit: db
    script:
    "checkv download_database checkv"
    stub:
    "mkdir -p checkv"
}

process DL_BAKTA {
    label 'process_medium'
    output:
    path 'bakta', emit: db
    script:
    "bakta_db download --output bakta --type full"
    stub:
    "mkdir -p bakta"
}

process DL_DRAM {
    label 'process_high'
    output:
    path 'dram', emit: db
    script:
    "DRAM-setup.py prepare_databases --output_dir dram --threads ${task.cpus}"
    stub:
    "mkdir -p dram"
}

process DL_SINGLEM {
    label 'process_low'
    output:
    path 'singlem', emit: db
    script:
    "singlem data --output-directory singlem"
    stub:
    "mkdir -p singlem"
}
