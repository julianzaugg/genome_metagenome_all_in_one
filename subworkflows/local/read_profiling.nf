/*
 * READ_PROFILING — sylph (+ sylph-tax) + singlem from raw reads.
 * sylph: per-sample sketch -> single combined profile -> taxonomy + merge.
 * singlem: one multi-sample pipe call over all samples.
 */

include { SYLPH_SKETCH; SYLPH_PROFILE; SYLPH_TAX } from '../../modules/local/sylph'
include { SINGLEM_PIPE }                           from '../../modules/local/singlem'

workflow READ_PROFILING {
    take:
    reads               // [ meta, reads ]
    sylph_db            // collected list of .syldb files
    sylph_tax_metadata  // collected list of *_metadata.tsv.gz (may be empty)
    singlem_metapackage // path
    run_sylph           // bool
    run_sylph_tax       // bool (taxonomy step; requires metadata)
    run_singlem         // bool

    main:
    ch_sylph_profile = Channel.empty()
    ch_sylph_rel     = Channel.empty()
    ch_sylph_seq     = Channel.empty()
    ch_singlem       = Channel.empty()
    ch_read_sets     = reads.collect(flat: false)

    if (run_sylph) {
        SYLPH_SKETCH(reads)
        ch_sketches = SYLPH_SKETCH.out.sketch.map { meta, s -> s }.collect()
        SYLPH_PROFILE(ch_sketches, sylph_db)
        ch_sylph_profile = SYLPH_PROFILE.out.profile

        if (run_sylph_tax) {
            SYLPH_TAX(SYLPH_PROFILE.out.profile, sylph_tax_metadata)
            ch_sylph_rel = SYLPH_TAX.out.relative
            ch_sylph_seq = SYLPH_TAX.out.sequence
        }
    }

    if (run_singlem) {
        ch_singlem_forward = ch_read_sets.map { rows ->
            rows.collect { row ->
                def meta = row[0]
                def r    = row[1]
                meta.single_end ? r : r[0]
            }
        }
        ch_singlem_reverse = ch_read_sets.map { rows ->
            rows.findAll { row -> !row[0].single_end }.collect { row -> row[1][1] }
        }
        SINGLEM_PIPE(ch_singlem_forward, ch_singlem_reverse, singlem_metapackage)
        ch_singlem = SINGLEM_PIPE.out.profile
    }

    emit:
    sylph              = ch_sylph_profile
    sylph_relative     = ch_sylph_rel
    sylph_sequence     = ch_sylph_seq
    singlem            = ch_singlem
}
