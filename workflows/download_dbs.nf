/*
 * DOWNLOAD_DBS — fetch/stage reference databases into params.db_outdir.
 *
 * Select with --download_db <name,name,...|all>. Users can instead point the
 * --*_db params at existing locations and never run this. See docs/databases.md.
 */

include { DL_GTDBTK; DL_CHECKM2; DL_GENOMAD; DL_CHECKV;
          DL_BAKTA;  DL_DRAM;    DL_SINGLEM } from '../modules/local/download_dbs'

workflow DOWNLOAD_DBS {

    def requested = (params.download_db ?: 'all').toString().toLowerCase()
                        .split(',').collect { it.trim() } as Set
    def want = { String name -> requested.contains('all') || requested.contains(name) }

    log.info "[gmaio] download_dbs: ${requested.join(', ')} -> ${params.outdir}/${params.db_outdir}"

    if (want('gtdbtk'))  DL_GTDBTK()
    if (want('checkm2')) DL_CHECKM2()
    if (want('genomad')) DL_GENOMAD()
    if (want('checkv'))  DL_CHECKV()
    if (want('bakta'))   DL_BAKTA()
    if (want('dram'))    DL_DRAM()
    if (want('singlem')) DL_SINGLEM()
}
