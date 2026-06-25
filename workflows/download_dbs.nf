/*
 * DOWNLOAD_DBS — SCAFFOLD (v1: wiring + TODO stubs)
 *
 * Fetch/stage reference databases into params.db_outdir using each tool's own
 * downloader, so users can either point the --*_db params at existing locations
 * OR build them here. Select with --download_db <name,name,...|all>.
 *
 *   gtdbtk   -> `gtdbtk download-db`
 *   genomad  -> `genomad download-database`
 *   checkm2  -> `checkm2 database --download`
 *   checkv   -> `checkv download_database`
 *   bakta    -> `bakta_db download`
 *   dram     -> `DRAM-setup.py prepare_databases`
 *   singlem  -> fetch metapackage
 */

workflow DOWNLOAD_DBS {

    def requested = (params.download_db ?: 'all').toString().toLowerCase()
    log.info "[gmaio] download_dbs: requested = '${requested}' -> ${params.db_outdir}"

    // TODO: one process per DB (label process_low, own container), guarded by `requested`.
    //       Each writes to "${params.db_outdir}/<name>" and prints the path to set as --<name>_db.

    log.warn "[gmaio] mode 'download_dbs' is a SCAFFOLD — DB download processes are not yet implemented. See docs/databases.md."
}
