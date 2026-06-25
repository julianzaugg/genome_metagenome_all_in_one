/*
 * Pool geNomad virus/plasmid sequences across samples, namespacing headers by
 * sample (derived from each file's geNomad prefix), ready for CheckV + clustering.
 * Mirrors process_genomad_scaffolds.sh.
 */

process GENOMAD_PROCESS {
    label 'process_low'

    input:
    path(virus_fnas,   stageAs: 'virus/*')    // <sample>_virus.fna.gz per sample
    path(plasmid_fnas, stageAs: 'plasmid/*')  // <sample>_plasmid.fna.gz per sample

    output:
    path 'genomad_virus.fna',   emit: virus
    path 'genomad_plasmid.fna', emit: plasmid
    path 'versions.yml',        emit: versions

    script:
    """
    namespace() {
        local indir=\$1 out=\$2 tag=\$3
        : > \$out
        for f in \$indir/*; do
            [ -e "\$f" ] || continue
            sample=\$(basename "\$f" | sed "s/_\${tag}.fna.*//")
            ( [[ "\$f" == *.gz ]] && zcat "\$f" || cat "\$f" ) \\
                | sed "s/^>/>\${sample}__genomad__/g" >> \$out
        done
    }
    namespace virus   genomad_virus.fna   virus
    namespace plasmid genomad_plasmid.fna plasmid

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo ">s__genomad__v1" > genomad_virus.fna;   echo "ACGT" >> genomad_virus.fna
    echo ">s__genomad__p1" > genomad_plasmid.fna; echo "ACGT" >> genomad_plasmid.fna
    echo '"${task.process}": {seqkit: stub}' > versions.yml
    """
}
