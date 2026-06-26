/*
 * Host / contaminant read removal with cleanifier (k-mer based).
 * Mirrors run_cleanifier.sh; paired-end uses --fastq R1 --pairs R2 per
 * https://gitlab.com/rahmannlab/cleanifier . Kept (host-removed) reads are
 * renamed to standard names for downstream steps.
 */

process CLEANIFIER_INDEX {
    tag   { host_ref.simpleName }
    label 'process_high'

    input:
    path(host_ref)
    val(nobjects)

    output:
    path 'cleanifier_host.*', emit: index
    path 'versions.yml',      emit: versions

    script:
    def args = task.ext.args ?: ''
    def nobjects_opt = nobjects ?: ''
    """
    if [ -n "${nobjects_opt}" ]; then
        nobjects="${nobjects_opt}"
    else
        nobjects=\$(python - "${host_ref}" <<'PY'
import bz2
import gzip
import lzma
import sys
from pathlib import Path

path = Path(sys.argv[1])
name = path.name.lower()
if name.endswith('.gz'):
    opener = gzip.open
elif name.endswith(('.bz2', '.bzip2')):
    opener = bz2.open
elif name.endswith(('.xz', '.lzma')):
    opener = lzma.open
elif name.endswith('.zst'):
    raise SystemExit('Cannot estimate --cleanifier_nobjects for .zst FASTA; provide --cleanifier_nobjects explicitly.')
else:
    opener = open

bases = 0
with opener(path, 'rt', errors='ignore') as handle:
    for line in handle:
        if not line.startswith('>'):
            bases += len(line.strip())

print(max(bases, 1))
PY
)
    fi

    cleanifier index ${args} \\
        --index cleanifier_host \\
        --files "${host_ref}" \\
        -n "\${nobjects}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cleanifier: \$(cleanifier --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch cleanifier_host.filter cleanifier_host.info
    echo '"${task.process}": {cleanifier: stub}' > versions.yml
    """
}

process CLEANIFIER {
    tag   { meta.id }
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta.id}.clean*.fastq.gz"), emit: reads
    path 'versions.yml',                                 emit: versions

    script:
    def args  = task.ext.args ?: ''
    def input = meta.single_end ? "--fastq ${reads}" : "--fastq ${reads[0]} --pairs ${reads[1]}"
    def index_file = index instanceof Collection ? index.find { it.name.endsWith('.filter') } : index
    if (!index_file) {
        throw new IllegalArgumentException("Cleanifier index input must include a .filter file")
    }
    """
    cleanifier filter ${args} \\
        --threads ${task.cpus} \\
        --index "${index_file}" \\
        ${input} \\
        --prefix ${meta.id}

    # cleanifier writes kept reads as <prefix>...keep...{,.1/.2}.f(ast)q.gz —
    # standardise names (adjust the glob if your cleanifier version differs).
    if [ "${meta.single_end}" = "true" ]; then
        mv \$(ls ${meta.id}*keep*.f*q.gz | head -1) ${meta.id}.clean.fastq.gz
    else
        keeps=(\$(ls ${meta.id}*keep*.f*q.gz | sort))
        mv "\${keeps[0]}" ${meta.id}.clean_1.fastq.gz
        mv "\${keeps[1]}" ${meta.id}.clean_2.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cleanifier: \$(cleanifier --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    if [ "${meta.single_end}" = "true" ]; then
        echo | gzip > ${meta.id}.clean.fastq.gz
    else
        echo | gzip > ${meta.id}.clean_1.fastq.gz
        echo | gzip > ${meta.id}.clean_2.fastq.gz
    fi
    echo '"${task.process}": {cleanifier: stub}' > versions.yml
    """
}
