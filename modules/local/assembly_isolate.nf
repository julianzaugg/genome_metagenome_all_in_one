/*
 * Isolate/metagenome assembly wrappers not provided by installed nf-core modules.
 */

process SHOVILL {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}"),                 emit: outdir
    path 'versions.yml',                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        --cpus ${task.cpus} \\
        --ram ${task.memory.toGiga()} \\
        --outdir ${meta.id} \\
        ${args}
    cp ${meta.id}/contigs.fa ${meta.id}.scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(shovill --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    echo ">${meta.id}_contig1" > ${meta.id}.scaffolds.fasta
    echo "ACGTACGT" >> ${meta.id}.scaffolds.fasta
    cp ${meta.id}.scaffolds.fasta ${meta.id}/contigs.fa
    echo '"${task.process}": {shovill: stub}' > versions.yml
    """
}

process MYLOASM {
    tag   { meta.id }
    label 'process_maxmem'
    label 'process_long'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}"),                 emit: outdir
    path 'versions.yml',                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${meta.id}
    myloasm ${args} \\
        --reads ${reads} \\
        --threads ${task.cpus} \\
        --output ${meta.id}

    assembly=\$(find ${meta.id} -type f \\( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \\) | head -1)
    if [[ -z "\$assembly" ]]; then
        echo "MYLOASM did not produce a FASTA assembly" >&2
        exit 1
    fi
    cp "\$assembly" ${meta.id}.scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        myloasm: \$(myloasm --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    echo ">${meta.id}_contig1" > ${meta.id}.scaffolds.fasta
    echo "ACGTACGT" >> ${meta.id}.scaffolds.fasta
    cp ${meta.id}.scaffolds.fasta ${meta.id}/assembly.fasta
    echo '"${task.process}": {myloasm: stub}' > versions.yml
    """
}

process AUTOCYCLER_ASSEMBLE {
    tag   { meta.id }
    label 'process_maxmem'
    label 'process_long'

    input:
    tuple val(meta), path(reads)
    val read_type

    output:
    tuple val(meta), path("${meta.id}.autocycler.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}"),                  emit: outdir
    path 'versions.yml',                                  emit: versions

    script:
    def args = task.ext.args ?: ''
    def jobs = task.ext.jobs ?: 4
    """
    mkdir -p ${meta.id}
    cd ${meta.id}
    genome_size=\$(autocycler helper genome_size --reads ../${reads} --threads ${task.cpus})
    autocycler subsample --reads ../${reads} --out_dir subsampled_reads --genome_size "\$genome_size"

    mkdir -p assemblies
    rm -f assemblies/jobs.txt
    for assembler in raven myloasm miniasm flye metamdbg necat nextdenovo plassembler canu; do
        for i in 01 02 03 04; do
            echo "autocycler helper \$assembler --reads subsampled_reads/sample_\$i.fastq --out_prefix assemblies/\${assembler}_\$i --threads ${task.cpus} --genome_size \$genome_size --read_type ${read_type}" >> assemblies/jobs.txt
        done
    done
    parallel --jobs ${jobs} --joblog assemblies/joblog.tsv --results assemblies/logs < assemblies/jobs.txt
    rm -f subsampled_reads/*.fastq

    autocycler compress -i assemblies -a autocycler_out ${args}
    autocycler cluster -a autocycler_out
    for c in autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c "\$c"
        autocycler resolve -c "\$c"
    done
    autocycler combine -a autocycler_out -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa
    cp autocycler_out/consensus_assembly.fasta ../${meta.id}.autocycler.fasta

    cd ..
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    echo ">${meta.id}_contig1" > ${meta.id}.autocycler.fasta
    echo "ACGTACGT" >> ${meta.id}.autocycler.fasta
    echo '"${task.process}": {autocycler: stub}' > versions.yml
    """
}

process POLYPOLISH {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(assembly), path(reads)

    output:
    tuple val(meta), path("${meta.id}.polypolished.fasta"), emit: assembly
    path 'versions.yml',                                    emit: versions

    when:
    meta.has_short_reads

    script:
    def args = task.ext.args ?: ''
    """
    minimap2 -ax sr -t ${task.cpus} ${assembly} ${reads[0]} ${reads[1]} > read_alignments.sam
    samtools view -@ ${task.cpus} -h -f 0x40 read_alignments.sam > alignments_1.sam
    samtools view -@ ${task.cpus} -h -f 0x80 read_alignments.sam > alignments_2.sam

    polypolish filter \\
        --in1 alignments_1.sam \\
        --in2 alignments_2.sam \\
        --out1 filtered_1.sam \\
        --out2 filtered_2.sam
    polypolish polish ${args} ${assembly} filtered_1.sam filtered_2.sam > ${meta.id}.polypolished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$(polypolish --version 2>&1 | head -1)
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    echo ">${meta.id}" > ${meta.id}.polypolished.fasta
    echo "ACGT" >> ${meta.id}.polypolished.fasta
    echo '"${task.process}": {polypolish: stub}' > versions.yml
    """
}

process DNAAPLER {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}"),                 emit: outdir
    path 'versions.yml',                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    dnaapler all \\
        --threads ${task.cpus} \\
        --prefix ${meta.id} \\
        --output ${meta.id} \\
        --input ${assembly} \\
        ${args}
    out=\$(find ${meta.id} -type f \\( -name '*reoriented*.fasta' -o -name '*.fasta' \\) | head -1)
    cp "\$out" ${meta.id}.scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: \$(dnaapler --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    echo ">${meta.id}" > ${meta.id}.scaffolds.fasta
    echo "ACGT" >> ${meta.id}.scaffolds.fasta
    cp ${meta.id}.scaffolds.fasta ${meta.id}/${meta.id}.fasta
    echo '"${task.process}": {dnaapler: stub}' > versions.yml
    """
}
