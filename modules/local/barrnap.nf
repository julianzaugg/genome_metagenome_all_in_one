/*
 * Barrnap — rRNA gene prediction (per genome/bin).
 * Mirrors run_barrnap.sh core: barrnap --outseq (gff + rRNA fasta); also pulls
 * out the 16S sequences for convenience. The elaborate per-hit length/coverage
 * table in the bash script is project-specific and omitted.
 */

process BARRNAP {
    tag   { meta.id }
    label 'process_low'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.rRNA.gff"),       emit: gff
    tuple val(meta), path("${meta.id}.rRNA.fasta"),     emit: rrna
    tuple val(meta), path("${meta.id}.16S.fasta"),      emit: ssu, optional: true
    path 'versions.yml',                                emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    barrnap ${args} --threads ${task.cpus} --quiet \\
        --outseq ${meta.id}.rRNA.fasta ${fasta} > ${meta.id}.rRNA.gff

    # extract 16S records (barrnap names them '>16S_rRNA::...')
    awk '/^>/{keep=(\$0 ~ /^>16S_rRNA/)} keep' ${meta.id}.rRNA.fasta > ${meta.id}.16S.fasta || true
    [ -s ${meta.id}.16S.fasta ] || rm -f ${meta.id}.16S.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        barrnap: \$(barrnap --version 2>&1 | sed 's/barrnap //')
    END_VERSIONS
    """

    stub:
    """
    echo "##gff-version 3" > ${meta.id}.rRNA.gff
    echo ">16S_rRNA::${meta.id}" > ${meta.id}.rRNA.fasta; echo "ACGT" >> ${meta.id}.rRNA.fasta
    cp ${meta.id}.rRNA.fasta ${meta.id}.16S.fasta
    echo '"${task.process}": {barrnap: stub}' > versions.yml
    """
}
