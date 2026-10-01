/*
 * Small utility processes.
 * PREP_ASSEMBLY — decompress + standardise an assembly to <id>.scaffolds.fasta
 * (nf-core SPADES emits *.scaffolds.fa.gz; downstream tools want plain fasta),
 * dropping short scaffolds (seqkit seq ext.args, e.g. '--min-len 500').
 */

process PREP_ASSEMBLY {
    tag   { meta.id }
    label 'process_single'

    input:
    tuple val(meta), path(scaffolds)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: assembly
    path 'versions.yml',                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    seqkit seq ${args} ${scaffolds} > ${meta.id}.scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.scaffolds.fasta
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}

/*
 * COLLECT_SINGLETONS — fastp's unpaired survivors (a read whose mate failed QC;
 * --unpaired1/--unpaired2, emitted by nf-core FASTP as *_R{1,2}.fail.fastq.gz)
 * concatenated into one file for metaSPAdes -s. gzip members concatenate cleanly.
 */
process COLLECT_SINGLETONS {
    tag   { meta.id }
    label 'process_single'

    input:
    tuple val(meta), path(fail_reads)

    output:
    tuple val(meta), path("${meta.id}.singletons.fastq.gz"), emit: reads

    script:
    """
    cat ${meta.id}_R1.fail.fastq.gz ${meta.id}_R2.fail.fastq.gz > ${meta.id}.singletons.fastq.gz
    """

    stub:
    """
    echo -e "@s1\\nACGT\\n+\\nIIII" | gzip > ${meta.id}.singletons.fastq.gz
    """
}
