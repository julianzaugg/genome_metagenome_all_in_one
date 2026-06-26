/*
 * HOST_REMOVAL — strip host/contaminant reads with cleanifier.
 * Handles both paired short reads (--fastq R1 --pairs R2) and single/long reads
 * (--fastq) per-sample via meta.single_end. Requires a cleanifier index
 * (--cleanifier_db).
 */

include { CLEANIFIER } from '../../modules/local/host_removal'

workflow HOST_REMOVAL {
    take:
    reads             // [ meta, reads ]
    cleanifier_index  // path to the cleanifier .filter index

    main:
    CLEANIFIER(reads, cleanifier_index)

    emit:
    reads = CLEANIFIER.out.reads
}
