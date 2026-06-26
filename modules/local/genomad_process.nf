/*
 * Pool geNomad virus/plasmid outputs across samples, namespacing by sample
 * (derived from each file's geNomad prefix). Mirrors process_genomad_scaffolds.sh:
 * pools sequences (.fna), proteins (.faa), gene tables and summary tables.
 * The .fna outputs feed CheckV + clustering; the rest are provenance records.
 */

process GENOMAD_PROCESS {
    label 'process_low'

    input:
    path(virus_fna,      stageAs: 'in/vfna/*')
    path(plasmid_fna,    stageAs: 'in/pfna/*')
    path(virus_faa,      stageAs: 'in/vfaa/*')
    path(plasmid_faa,    stageAs: 'in/pfaa/*')
    path(virus_genes,    stageAs: 'in/vgenes/*')
    path(plasmid_genes,  stageAs: 'in/pgenes/*')
    path(virus_summary,  stageAs: 'in/vsum/*')
    path(plasmid_summary,stageAs: 'in/psum/*')

    output:
    path 'genomad_virus.fna',            emit: virus
    path 'genomad_plasmid.fna',          emit: plasmid
    path 'genomad_virus_proteins.faa',   emit: virus_proteins,   optional: true
    path 'genomad_plasmid_proteins.faa', emit: plasmid_proteins, optional: true
    path 'genomad_virus_genes.tsv',      emit: virus_genes,      optional: true
    path 'genomad_plasmid_genes.tsv',    emit: plasmid_genes,    optional: true
    path 'genomad_virus_summary.tsv',    emit: virus_summary,    optional: true
    path 'genomad_plasmid_summary.tsv',  emit: plasmid_summary,  optional: true
    path 'versions.yml',                 emit: versions

    script:
    """
    # \$1 dir, \$2 out, \$3 suffix-to-strip — concatenate FASTA, prefix headers by sample
    namespace_seq() {
        : > "\$2"
        for f in "\$1"/*; do
            [ -e "\$f" ] || continue
            s=\$(basename "\$f" | sed "s/\$3.*//")
            ( [[ "\$f" == *.gz ]] && zcat "\$f" || cat "\$f" ) | sed "s/^>/>\${s}__genomad__/g" >> "\$2"
        done
    }
    # \$1 dir, \$2 out, \$3 suffix — header from first file, prefix data rows by sample
    namespace_tsv() {
        : > "\$2"; local hdr=
        for f in "\$1"/*; do
            [ -e "\$f" ] || continue
            s=\$(basename "\$f" | sed "s/\$3.*//")
            [ -z "\$hdr" ] && { head -n 1 "\$f" > "\$2"; hdr=1; }
            tail -n +2 "\$f" | sed "s/^/\${s}__genomad__/g" >> "\$2"
        done
    }

    namespace_seq in/vfna   genomad_virus.fna            _virus.fna
    namespace_seq in/pfna   genomad_plasmid.fna          _plasmid.fna
    namespace_seq in/vfaa   genomad_virus_proteins.faa   _virus_proteins.faa
    namespace_seq in/pfaa   genomad_plasmid_proteins.faa _plasmid_proteins.faa
    namespace_tsv in/vgenes genomad_virus_genes.tsv      _virus_genes.tsv
    namespace_tsv in/pgenes genomad_plasmid_genes.tsv    _plasmid_genes.tsv
    namespace_tsv in/vsum   genomad_virus_summary.tsv    _virus_summary.tsv
    namespace_tsv in/psum   genomad_plasmid_summary.tsv  _plasmid_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo ">s__genomad__v1" > genomad_virus.fna;   echo "ACGT" >> genomad_virus.fna
    echo ">s__genomad__p1" > genomad_plasmid.fna; echo "ACGT" >> genomad_plasmid.fna
    touch genomad_virus_proteins.faa genomad_plasmid_proteins.faa
    touch genomad_virus_genes.tsv genomad_plasmid_genes.tsv
    touch genomad_virus_summary.tsv genomad_plasmid_summary.tsv
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}
