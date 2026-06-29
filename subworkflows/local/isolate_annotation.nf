/*
 * ISOLATE_ANNOTATION — Bakta, MLST, AMRFinderPlus and ISEScan for isolate
 * assemblies.
 */

include {
    BAKTA_BAKTA;
    BAKTA_STATS;
    MLST;
    AMRFINDERPLUS_RUN;
    AMRFINDERPLUS_COLLATE;
    ISESCAN;
    ISESCAN_COLLATE
} from '../../modules/local/isolate_annotation'

workflow ISOLATE_ANNOTATION {
    take:
    assemblies               // [ meta, fasta ]
    bakta_db                 // path
    bakta_reference_proteins // path or []
    amrfinder_db             // path or []
    run_annotation           // bool
    run_amrfinder            // bool
    run_isescan              // bool

    main:
    ch_gff     = Channel.empty()
    ch_faa     = Channel.empty()
    ch_json    = Channel.empty()
    ch_stats   = Channel.empty()
    ch_mlst    = Channel.empty()
    ch_amr     = Channel.empty()
    ch_isescan = Channel.empty()
    ch_versions = Channel.empty()

    if (run_annotation) {
        BAKTA_BAKTA(assemblies, bakta_db, bakta_reference_proteins)
        ch_gff = BAKTA_BAKTA.out.gff
        ch_faa = BAKTA_BAKTA.out.faa
        ch_json = BAKTA_BAKTA.out.json
        BAKTA_STATS(ch_json.map { meta, json -> json }.collect())
        ch_stats = BAKTA_STATS.out.stats
        ch_versions = ch_versions
            .mix(BAKTA_BAKTA.out.versions)
            .mix(BAKTA_STATS.out.versions)

        MLST(assemblies.map { meta, fasta -> fasta }.collect())
        ch_mlst = MLST.out.summary
        ch_versions = ch_versions.mix(MLST.out.versions)

        if (run_amrfinder) {
            ch_amr_in = ch_faa.join(ch_gff)
            AMRFINDERPLUS_RUN(ch_amr_in, amrfinder_db)
            AMRFINDERPLUS_COLLATE(AMRFINDERPLUS_RUN.out.results.map { meta, t -> t }.collect())
            ch_amr = AMRFINDERPLUS_COLLATE.out.summary
            ch_versions = ch_versions
                .mix(AMRFINDERPLUS_RUN.out.versions)
                .mix(AMRFINDERPLUS_COLLATE.out.versions)
        }
    }

    if (run_isescan) {
        ISESCAN(assemblies)
        ISESCAN_COLLATE(ISESCAN.out.tables.map { meta, t -> t }.collect())
        ch_isescan = ISESCAN_COLLATE.out.summary
        ch_versions = ch_versions
            .mix(ISESCAN.out.versions)
            .mix(ISESCAN_COLLATE.out.versions)
    }

    emit:
    gff       = ch_gff
    faa       = ch_faa
    json      = ch_json
    stats     = ch_stats
    mlst      = ch_mlst
    amrfinder = ch_amr
    isescan   = ch_isescan
    versions  = ch_versions
}
