/*
 * RPKM: length-filter selected Illumina R1 reads, DIAMOND blastx against the
 * gene catalogue and SingleM marker databases, then calculate normalized RPKM.
 */

process RPKM_FILTER_READS {
    tag   { meta.id }
    label 'process_low'

    input:
    tuple val(meta), path(read)
    val(min_len)

    output:
    tuple val(meta), path("${meta.id}.rpkm.min${min_len}.fastq.gz"), emit: reads
    path 'versions.yml', emit: versions

    script:
    """
    seqkit seq --min-len ${min_len} --threads ${task.cpus} \\
        --out-file ${meta.id}.rpkm.min${min_len}.fastq.gz \\
        ${read}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/^seqkit //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.rpkm.min${min_len}.fastq.gz
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}

process RPKM_SINGLEM_MARKERS {
    tag   { metapackage.simpleName }
    label 'process_low'

    input:
    path(metapackage)

    output:
    path 'singlem_marker_diamond',       emit: marker_dbs
    path 'singlem_marker_lengths.tsv',   emit: marker_lengths
    path 'versions.yml',                 emit: versions

    script:
    """
    set -euo pipefail

    mkdir -p singlem_marker_diamond

    package_root="${metapackage}"
    if [ -d "${metapackage}/payload_directory" ]; then
        package_root="${metapackage}/payload_directory"
    fi

    found_dmnd=0
    for spkg in "\${package_root}"/*.spkg; do
        [ -d "\${spkg}" ] || continue
        package_name=\$(basename "\${spkg}" .spkg)
        dmnd="\${spkg}/\${package_name}/S3.dmnd"
        if [ -f "\${dmnd}" ]; then
            ln -s "\${dmnd}" "singlem_marker_diamond/\${package_name}.dmnd"
            found_dmnd=1
        fi
    done

    if [ "\${found_dmnd}" -eq 0 ]; then
        echo "[RPKM_SINGLEM_MARKERS] No SingleM marker S3.dmnd files found under \${package_root}" >&2
        exit 1
    fi

    echo -e "Gene\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len" > singlem_marker_lengths.tsv
    found_faa=0
    while IFS= read -r faa; do
        found_faa=1
        spkg=\$(printf '%s\\n' "\${faa}" | awk -F/ '{ for (i=1; i<=NF; i++) if (\$i ~ "\\\\.spkg\$") print \$i }' | tail -n 1)
        gene=\${spkg%.spkg}
        seqkit stats --tabular "\${faa}" \\
            | awk -v gene="\${gene}" 'BEGIN{OFS="\\t"} NR==2{print gene,\$4,\$5,\$6,\$7,\$8}' \\
            >> singlem_marker_lengths.tsv
    done < <(find "\${package_root}" -type f -path "*.spkg/*/*.faa" | sort)
    if [ "\${found_faa}" -eq 0 ]; then
        echo "[RPKM_SINGLEM_MARKERS] No SingleM marker FASTA files found under \${package_root}" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/^seqkit //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p singlem_marker_diamond
    touch singlem_marker_diamond/S3.stub.dmnd
    echo -e "Gene\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len" > singlem_marker_lengths.tsv
    echo -e "S3.stub\\t1\\t300\\t100\\t100\\t100" >> singlem_marker_lengths.tsv
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}

process RPKM_DIAMOND_MAKEDB {
    label 'process_low'

    input:
    path(catalogue)

    output:
    path 'gene_catalogue.dmnd', emit: db
    path 'versions.yml',        emit: versions

    script:
    """
    diamond makedb --in ${catalogue} --db gene_catalogue --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch gene_catalogue.dmnd
    echo '"${task.process}": {diamond: stub}' > versions.yml
    """
}

process RPKM_DIAMOND_GENE {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(read), path(db)

    output:
    tuple val(meta), path("${meta.id}.gene_catalogue_blast.tsv"), emit: blast
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: '--evalue 0.00001 --min-score 40 --query-cover 80 --id 70 --max-hsps 1 --max-target-seqs 1'
    """
    tmp=${meta.id}.gene_catalogue_blast.tmp
    diamond blastx \\
        --query ${read} \\
        --db ${db} \\
        ${args} \\
        --threads ${task.cpus} \\
        --outfmt 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \\
        --out "\${tmp}"

    echo -e "sample\\tqseqid\\tsseqid\\tstitle\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tpercent_query_aligned\\tpercent_subject_aligned" > ${meta.id}.gene_catalogue_blast.tsv
    awk -v sample="${meta.id}" -v OFS='\\t' -F '\\t' 'function abs(v){return v < 0 ? -v : v} {pq=(\$14 != 0) ? sprintf("%.3f", abs(\$9-\$8)/\$14) : "NA"; ps=(\$15 != 0) ? sprintf("%.3f", abs(\$11-\$10)/\$15) : "NA"; print sample,\$0,pq,ps}' "\${tmp}" >> ${meta.id}.gene_catalogue_blast.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    echo -e "sample\\tqseqid\\tsseqid\\tstitle\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tpercent_query_aligned\\tpercent_subject_aligned" > ${meta.id}.gene_catalogue_blast.tsv
    echo -e "${meta.id}\\tread1\\tgene1\\tgene1\\t100\\t100\\t0\\t0\\t1\\t100\\t1\\t100\\t1e-20\\t100\\t150\\t300\\t0.660\\t0.330" >> ${meta.id}.gene_catalogue_blast.tsv
    echo '"${task.process}": {diamond: stub}' > versions.yml
    """
}

process RPKM_COLLATE_GENE_BLAST {
    label 'process_single'

    input:
    path(blasts)

    output:
    path 'gene_catalogue_blast.tsv', emit: blast
    path 'versions.yml',             emit: versions

    script:
    """
    first=1
    for blast in ${blasts}; do
        if [ "\${first}" -eq 1 ]; then
            cat "\${blast}" > gene_catalogue_blast.tsv
            first=0
        else
            tail -n +2 "\${blast}" >> gene_catalogue_blast.tsv
        fi
    done

    echo '"${task.process}": {bash: true}' > versions.yml
    """

    stub:
    """
    echo -e "sample\\tqseqid\\tsseqid\\tstitle\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tpercent_query_aligned\\tpercent_subject_aligned" > gene_catalogue_blast.tsv
    echo -e "SAMPLE_A\\tread1\\tgene1\\tgene1\\t100\\t100\\t0\\t0\\t1\\t100\\t1\\t100\\t1e-20\\t100\\t150\\t300\\t0.660\\t0.330" >> gene_catalogue_blast.tsv
    echo '"${task.process}": {bash: stub}' > versions.yml
    """
}

process RPKM_DIAMOND_SINGLEM {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(read), path(marker_dbs)

    output:
    path "${meta.id}.singlem_marker_blast", emit: blast_dir
    path 'versions.yml', emit: versions

    script:
    def args = task.ext.args ?: '--evalue 0.00001 --min-score 40 --query-cover 80 --id 70 --max-hsps 1 --max-target-seqs 1'
    """
    mkdir -p ${meta.id}.singlem_marker_blast
    for db in ${marker_dbs}/*.dmnd; do
        marker=\$(basename "\${db}" .dmnd)
        tmp="${meta.id}.singlem_marker_blast/\${marker}.tmp"
        out="${meta.id}.singlem_marker_blast/\${marker}_blast.tsv"

        diamond blastx \\
            --query ${read} \\
            --db "\${db}" \\
            ${args} \\
            --threads ${task.cpus} \\
            --outfmt 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \\
            --out "\${tmp}"

        echo -e "sample\\tqseqid\\tsseqid\\tstitle\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tpercent_query_aligned\\tpercent_subject_aligned" > "\${out}"
        awk -v sample="${meta.id}" -v OFS='\\t' -F '\\t' 'function abs(v){return v < 0 ? -v : v} {pq=(\$14 != 0) ? sprintf("%.3f", abs(\$9-\$8)/\$14) : "NA"; ps=(\$15 != 0) ? sprintf("%.3f", abs(\$11-\$10)/\$15) : "NA"; print sample,\$0,pq,ps}' "\${tmp}" >> "\${out}"
        rm -f "\${tmp}"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}.singlem_marker_blast
    echo -e "sample\\tqseqid\\tsseqid\\tstitle\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tpercent_query_aligned\\tpercent_subject_aligned" > ${meta.id}.singlem_marker_blast/S3.stub_blast.tsv
    echo -e "${meta.id}\\tread1\\tmarker1\\tmarker1\\t100\\t100\\t0\\t0\\t1\\t100\\t1\\t100\\t1e-20\\t100\\t150\\t300\\t0.660\\t0.330" >> ${meta.id}.singlem_marker_blast/S3.stub_blast.tsv
    echo '"${task.process}": {diamond: stub}' > versions.yml
    """
}

process RPKM_CALCULATE {
    label 'process_single'

    input:
    path(singlem_blast_dirs, stageAs: 'singlem_marker_blast/*')
    path(marker_lengths)
    path(gene_blast)

    output:
    path 'singlem_sample_rpkm.tsv',                     emit: singlem_rpkm
    path 'singlem_rpkm_means.tsv',                      emit: singlem_means
    path 'gene_catalogue_rpkm_per_gene_normalised.tsv', emit: normalised_rpkm
    path 'gene_catalogue_mapped_reads_per_gene.tsv',    emit: mapped_reads
    path 'versions.yml',                                emit: versions

    script:
    """
    python ${projectDir}/bin/rpkm_calculate.py \\
        --singlem-blast-dir singlem_marker_blast \\
        --marker-lengths ${marker_lengths} \\
        --gene-blast ${gene_blast} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "sample\\tMarker_gene\\trpkm" > singlem_sample_rpkm.tsv
    echo -e "SAMPLE_A\\tS3.stub\\t3333333.333333" >> singlem_sample_rpkm.tsv
    echo -e "sample\\tMean_rpkm" > singlem_rpkm_means.tsv
    echo -e "SAMPLE_A\\t3333333.333333" >> singlem_rpkm_means.tsv
    echo -e "Gene_ID\\tSAMPLE_A" > gene_catalogue_rpkm_per_gene_normalised.tsv
    echo -e "gene1\\t100" >> gene_catalogue_rpkm_per_gene_normalised.tsv
    echo -e "Gene_ID\\tSAMPLE_A" > gene_catalogue_mapped_reads_per_gene.tsv
    echo -e "gene1\\t1" >> gene_catalogue_mapped_reads_per_gene.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
