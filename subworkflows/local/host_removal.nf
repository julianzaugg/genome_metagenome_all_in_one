/*
 * HOST_REMOVAL — strip host/contaminant reads with cleanifier.
 * Handles both paired short reads (--fastq R1 --pairs R2) and single/long reads
 * (--fastq) per-sample via meta.single_end. Consumes a cleanifier index from
 * --cleanifier_db or one built from --host_ref.
 */

include { CLEANIFIER } from '../../modules/local/host_removal'

workflow HOST_REMOVAL {
    take:
    reads             // [ meta, reads ]
    cleanifier_index  // path to the cleanifier .filter index

    main:
    CLEANIFIER(reads, cleanifier_index)

    emit:
    reads    = CLEANIFIER.out.reads
    versions = CLEANIFIER.out.versions
}
