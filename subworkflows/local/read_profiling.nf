/*
 * READ_PROFILING — sylph + singlem from QC'd reads.
 * sylph: per-sample sketch then a single combined profile across all samples.
 * singlem: per-sample pipe.
 */

include { SYLPH_SKETCH; SYLPH_PROFILE } from '../../modules/local/sylph'
include { SINGLEM_PIPE }                from '../../modules/local/singlem'

workflow READ_PROFILING {
    take:
    reads               // [ meta, reads ]
    sylph_db            // path
    singlem_metapackage // path
    run_sylph           // bool
    run_singlem         // bool

    main:
    ch_sylph_profile  = Channel.empty()
    ch_singlem        = Channel.empty()

    if (run_sylph) {
        SYLPH_SKETCH(reads)
        ch_sketches = SYLPH_SKETCH.out.sketch.map { meta, s -> s }.collect()
        SYLPH_PROFILE(ch_sketches, sylph_db)
        ch_sylph_profile = SYLPH_PROFILE.out.profile
    }

    if (run_singlem) {
        SINGLEM_PIPE(reads, singlem_metapackage)
        ch_singlem = SINGLEM_PIPE.out.profile
    }

    emit:
    sylph   = ch_sylph_profile
    singlem = ch_singlem
}
