/*
 * COMPARISON_ASSEMBLY_PREP — decompress a comparison assembly if gzipped, and
 * rename it to the samplesheet's sample id (pyrodigal itself doesn't read
 * .gz directly; mirrors the decompression REFERENCE_PREP does for
 * --reference_genomes, modules/local/reference_genomes.nf).
 *
 * The input is staged under input/ so an assembly already named
 * <sample>.fasta can't collide with the output of the same name.
 */
process COMPARISON_ASSEMBLY_PREP {
    tag   { meta.id }
    label 'process_single'

    input:
    tuple val(meta), path(assembly, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: assembly
    path 'versions.yml',                       emit: versions

    script:
    """
    case "${assembly}" in
        *.gz) gzip -dc ${assembly} > ${meta.id}.fasta ;;
        *)    cp -L ${assembly} ${meta.id}.fasta ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //;s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo ">${meta.id}_1" > ${meta.id}.fasta
    echo "ACGT"          >> ${meta.id}.fasta
    echo '"${task.process}": {bash: stub}' > versions.yml
    """
}
