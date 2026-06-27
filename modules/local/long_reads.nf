/*
 * Long-read preprocessing and Dorado helpers.
 */

process DORADO_BASECALL {
    tag   { meta.id }
    label 'gpu'

    input:
    tuple val(meta), path(pod5_dir)
    val model
    val barcode_kit
    val device

    output:
    tuple val(meta), path("${meta.id}.basecalled.bam"), emit: bam
    path 'versions.yml',                               emit: versions

    script:
    def kit = barcode_kit ? "--kit-name ${barcode_kit}" : ''
    """
    dorado basecaller \\
        --device ${device} \\
        --recursive \\
        --min-qscore 0 \\
        --trim all \\
        ${kit} \\
        ${model} \\
        ${pod5_dir} > ${meta.id}.basecalled.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.basecalled.bam
    echo '"${task.process}": {dorado: stub}' > versions.yml
    """
}

process DORADO_DEMUX {
    tag   { meta.id }
    label 'gpu'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads
    path 'versions.yml',                          emit: versions

    script:
    """
    mkdir -p demux
    dorado demux \\
        --no-classify \\
        --emit-fastq \\
        --output-dir demux \\
        --threads ${task.cpus} \\
        ${bam}
    cat demux/*.fastq | gzip -n > ${meta.id}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo "@${meta.id}" | gzip -n > ${meta.id}.fastq.gz
    echo '"${task.process}": {dorado: stub}' > versions.yml
    """
}

process PORECHOP {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.porechop.fastq.gz"), emit: reads
    path 'versions.yml',                                  emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    porechop ${args} -i ${reads} -o ${meta.id}.porechop.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$(porechop --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.porechop.fastq.gz
    echo '"${task.process}": {porechop: stub}' > versions.yml
    """
}

process FASTPLONG {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.fastplong.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastplong.html"),     emit: html
    tuple val(meta), path("${meta.id}.fastplong.json"),     emit: json
    path 'versions.yml',                                    emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    fastplong \\
        -i ${reads} \\
        -o ${meta.id}.fastplong.fastq.gz \\
        -h ${meta.id}.fastplong.html \\
        -j ${meta.id}.fastplong.json \\
        --thread ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.fastplong.fastq.gz ${meta.id}.fastplong.html ${meta.id}.fastplong.json
    echo '"${task.process}": {fastplong: stub}' > versions.yml
    """
}

process DORADO_POLISH {
    tag   { meta.id }
    label 'gpu'

    input:
    tuple val(meta), path(assembly), path(reads)
    val device

    output:
    tuple val(meta), path("${meta.id}.dorado_polished.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}.dorado_variants.vcf"),   emit: variants, optional: true
    tuple val(meta), path("${meta.id}.dorado_aligned.bam"),    emit: bam,      optional: true
    path 'versions.yml',                                       emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    dorado aligner --threads ${task.cpus} ${assembly} ${reads} | \\
        samtools sort --threads ${task.cpus} > ${meta.id}.dorado_aligned.bam
    samtools index ${meta.id}.dorado_aligned.bam

    dorado polish ${args} \\
        ${meta.id}.dorado_aligned.bam \\
        ${assembly} \\
        --device ${device} \\
        --threads ${task.cpus} \\
        --infer-threads ${task.cpus} \\
        --gvcf \\
        -o dorado_polish

    mv dorado_polish/consensus.fasta ${meta.id}.dorado_polished.fasta
    if [[ -f dorado_polish/variants.vcf ]]; then
        mv dorado_polish/variants.vcf ${meta.id}.dorado_variants.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dorado: \$(dorado --version 2>&1 | head -1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    echo ">${meta.id}" > ${meta.id}.dorado_polished.fasta
    echo "ACGT" >> ${meta.id}.dorado_polished.fasta
    touch ${meta.id}.dorado_variants.vcf ${meta.id}.dorado_aligned.bam
    echo '"${task.process}": {dorado: stub}' > versions.yml
    """
}
