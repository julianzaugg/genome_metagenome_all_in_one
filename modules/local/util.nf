/*
 * Small utility processes.
 * PREP_ASSEMBLY — decompress + standardise an assembly to <id>.scaffolds.fasta
 * (nf-core SPADES emits *.scaffolds.fa.gz; downstream tools want plain fasta).
 */

process PREP_ASSEMBLY {
    tag   { meta.id }
    label 'process_single'

    input:
    tuple val(meta), path(scaffolds)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: assembly

    script:
    """
    if [[ "${scaffolds}" == *.gz ]]; then
        gzip -cd ${scaffolds} > ${meta.id}.scaffolds.fasta
    else
        cp ${scaffolds} ${meta.id}.scaffolds.fasta
    fi
    """

    stub:
    """
    touch ${meta.id}.scaffolds.fasta
    """
}
