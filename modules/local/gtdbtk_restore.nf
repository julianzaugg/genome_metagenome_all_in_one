/*
 * GTDBTK_RESTORE_NAMES — strip the reserved 'USERREF_' prefix from a GTDB-Tk
 * classify_wf output directory.
 *
 * External reference genomes are renamed 'USERREF_<stem>' before GTDB-Tk (so an
 * accession-style filename such as GCF_000123.1 cannot collide with a real GTDB
 * reference id and abort the run). This step produces a cleaned copy of the whole
 * output directory — summaries, classify trees and MSAs — with the original names
 * restored, so the published GTDB-Tk outputs read naturally. It is a no-op when no
 * references were classified (nothing carries the prefix).
 *
 * The marker-gene tree still reads the RAW (prefixed) directory and strips the
 * prefix itself when writing leaf labels, so this process only shapes what is
 * published.
 */

process GTDBTK_RESTORE_NAMES {
    tag   { meta.id }
    label 'process_single'

    input:
    tuple val(meta), path(gtdbtk_dir, stageAs: 'gtdbtk_in')

    output:
    tuple val(meta), path("${prefix}"),                           emit: gtdb_outdir
    tuple val(meta), path("${prefix}/classify/*.summary.tsv"),    emit: summary
    path 'versions.yml',                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -rL gtdbtk_in ${prefix}

    # Text outputs: summaries, classify trees, marker summaries. USERREF_ only ever
    # appears as a genome-id prefix, so a plain global substitution is safe.
    find ${prefix} -type f \\
        \\( -name '*.summary.tsv' -o -name '*.classify.tree' -o -name '*.markers_summary.tsv' \\) \\
        -exec sed -i.bak 's/USERREF_//g' {} \\;
    find ${prefix} -type f -name '*.bak' -delete

    # Gzipped MSAs (user_msa / msa): decompress, substitute, recompress in place.
    for gz in \$(find ${prefix} -type f -name '*.fasta.gz'); do
        gunzip -c "\$gz" | sed 's/USERREF_//g' | gzip > "\${gz}.tmp" && mv -f "\${gz}.tmp" "\$gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | grep -Eo '[0-9]+(\\.[0-9]+)+' | head -1 || echo NA)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -rL gtdbtk_in ${prefix}
    find ${prefix} -type f \\( -name '*.summary.tsv' -o -name '*.classify.tree' \\) \\
        -exec sed -i.bak 's/USERREF_//g' {} \\; 2>/dev/null || true
    find ${prefix} -type f -name '*.bak' -delete 2>/dev/null || true
    echo '"${task.process}": {sed: stub}' > versions.yml
    """
}
