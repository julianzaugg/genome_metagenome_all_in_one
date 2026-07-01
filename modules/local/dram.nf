/*
 * DRAM — functional annotation of the gene catalogue proteins.
 * Mirrors run_dram_cdhit.sh (DRAM.py annotate_genes on a protein fasta).
 */

process DRAM_ANNOTATE {
    label 'process_high'
    label 'process_long'

    input:
    path(proteins)
    path(dram_db)

    output:
    path 'dram_annotations/annotations.tsv', emit: annotations
    path 'dram_annotations',                 emit: outdir
    path 'versions.yml',                     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # DRAM reads DB locations from its config. Some installs keep CONFIG
    # outside the database directory, so only override it for prepared DB dirs
    # that actually carry their own CONFIG file.
    if [[ -f "${dram_db}/CONFIG" ]]; then
        export DRAM_CONFIG_LOCATION="${dram_db}/CONFIG"
    fi

    DRAM.py annotate_genes ${args} \\
        --input_faa ${proteins} \\
        --output_dir dram_annotations \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: \$(pip show dram 2>/dev/null | grep '^Version:' | sed 's/Version: //' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dram_annotations
    echo -e "\\tgene_id\\tannotation" > dram_annotations/annotations.tsv
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}

/*
 * DRAM — per-bin functional annotation. Mirrors run_dram_bins_parallel.sh
 * (DRAM.py annotate_genes on each bin's .faa). Nextflow's natural per-item
 * parallelism replaces GNU Parallel.
 */

process DRAM_ANNOTATE_BINS {
    tag { meta.id }
    label 'process_high'
    label 'process_long'

    input:
    tuple val(meta), path(proteins)
    path(dram_db)

    output:
    tuple val(meta), path('dram_annotations/annotations.tsv'), emit: annotations
    tuple val(meta), path('dram_annotations'),                  emit: outdir
    path 'versions.yml',                                        emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    if [[ -f "${dram_db}/CONFIG" ]]; then
        export DRAM_CONFIG_LOCATION="${dram_db}/CONFIG"
    fi

    DRAM.py annotate_genes ${args} \\
        --input_faa ${proteins} \\
        --output_dir dram_annotations \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: \$(pip show dram 2>/dev/null | grep '^Version:' | sed 's/Version: //' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dram_annotations
    echo -e "\\tgene_id\\tannotation" > dram_annotations/annotations.tsv
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}

process DRAM_DISTILL {
    tag { meta.id }
    label 'process_low'

    input:
    tuple val(meta), path(annotations)
    path(dram_db)

    output:
    tuple val(meta), path('distilled'), emit: distilled
    path 'versions.yml',                emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Build a local DRAM config with valid distill-sheet paths. The staged DB
    # CONFIG often has dram_sheets entries set to None (paths were absolute to
    # the setup host and don't resolve inside this container). annotate doesn't
    # need the sheets; distill does. The sheets ship with the mag_annotator
    # package (version-suffixed, e.g. genome_summary_form.<date>.tsv), so we
    # locate them by prefix and backfill any missing/None entries.
    python3 - <<'PYEOF'
import json, os, sys

# distill sheet keys DRAM needs (1.4.x)
KEYS = ['genome_summary_form', 'module_step_form', 'function_heatmap_form',
        'amg_database', 'etc_module_database']

src = '${dram_db}/CONFIG'
cfg = {}
if os.path.isfile(src):
    try:
        with open(src) as fh:
            cfg = json.load(fh)
    except Exception as e:
        sys.stderr.write("WARN: could not parse %s: %s\\n" % (src, e))
if not isinstance(cfg, dict):
    cfg = {}

try:
    import mag_annotator
    pkg = os.path.dirname(os.path.abspath(mag_annotator.__file__))
except Exception:
    pkg = ''

# Recursively collect candidate .tsv sheets from the package and staged DB dir.
found = {}
for root in [pkg, '${dram_db}']:
    if not root or not os.path.isdir(root):
        continue
    for dirpath, _dirs, files in os.walk(root):
        for fn in files:
            if not fn.endswith('.tsv'):
                continue
            for k in KEYS:
                if fn.startswith(k) and k not in found:
                    found[k] = os.path.join(dirpath, fn)

sh = cfg.setdefault('dram_sheets', {})
for k in KEYS:
    cur = sh.get(k)
    if (not cur or not os.path.isfile(str(cur))) and k in found:
        sh[k] = found[k]

sys.stderr.write("DRAM distill sheets resolved:\\n")
for k in KEYS:
    sys.stderr.write("  %s -> %s\\n" % (k, sh.get(k)))

with open('LOCAL_DRAM_CONFIG.json', 'w') as fh:
    json.dump(cfg, fh, indent=2)
PYEOF

    export DRAM_CONFIG_LOCATION="LOCAL_DRAM_CONFIG.json"

    DRAM.py distill ${args} \\
        --input_file ${annotations} \\
        --output_dir distilled

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dram: \$(pip show dram 2>/dev/null | grep '^Version:' | sed 's/Version: //' || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p distilled
    touch distilled/metabolism_summary.xlsx
    echo '"${task.process}": {dram: stub}' > versions.yml
    """
}
