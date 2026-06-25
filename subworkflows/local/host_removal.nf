/*
 * HOST_REMOVAL — strip host/contaminant reads.
 * Paired short reads -> MINIMAP2_HOSTFILTER (minimap2 -ax sr | samtools -f 13).
 * Single/long reads  -> CLEANIFIER.
 * Selection is per-sample on meta.single_end so a mixed run works.
 *
 * host_ref / cleanifier_index are resolved from params (one default each); a
 * per-sample host_ref column could be threaded here later.
 */

include { MINIMAP2_HOSTFILTER } from '../../modules/local/host_removal'
include { CLEANIFIER }          from '../../modules/local/host_removal'

workflow HOST_REMOVAL {
    take:
    reads             // [ meta, reads ]
    host_ref          // path (may be []),  used for paired short reads
    cleanifier_index  // path (may be []),  used for single/long reads

    main:
    branched = reads.branch {
        paired: !it[0].single_end
        single:  it[0].single_end
    }

    MINIMAP2_HOSTFILTER(branched.paired, host_ref)
    CLEANIFIER(branched.single, cleanifier_index)

    cleaned = MINIMAP2_HOSTFILTER.out.reads.mix(CLEANIFIER.out.reads)

    emit:
    reads = cleaned
}
