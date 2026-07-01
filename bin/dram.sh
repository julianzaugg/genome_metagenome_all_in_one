#!/usr/bin/env bash
#
# Wrapper around DRAM.py applying two runtime patches for DRAM container bugs we
# cannot fix in the read-only image. Config is supplied via DRAM_CONFIG_LOCATION.
#
# (1) DRAM 1.5.0 CLI bug: the annotate / annotate_genes / distill subparsers add
#     --config_loc (default None), but __main__ calls the handler with every
#     parsed arg except 'func', so config_loc=None is passed to a handler that
#     doesn't accept it and every invocation dies with TypeError. We run a local
#     copy of DRAM.py with config_loc also stripped from that dict. The sed is a
#     no-op on 1.4.6 (no such arg), so this is safe across versions.
#
# (2) DRAM 1.4.x annotate bug: generate_annotated_fasta accesses
#     annotation.pfam_hits unconditionally for genes graded rank "D", which
#     AttributeErrors when pfam was not run. We exclude pfam (its mmseqs profile
#     search crashes this container's mmseqs — see build_dram_config.py), so we
#     inject a sitecustomize.py that guards that access the same way the rank
#     "C" branch already guards its own columns. Harmless if pfam_hits exists.
set -euo pipefail

# (1) local, config_loc-stripped copy of DRAM.py
src=$(command -v DRAM.py)
cp "$src" ./.dram_local.py
sed -i "s/if i != 'func'/if i not in ('func', 'config_loc')/" ./.dram_local.py

# (2) sitecustomize patch, auto-imported at interpreter startup via PYTHONPATH
cat > sitecustomize.py <<'PYEOF'
try:
    import inspect
    import mag_annotator.annotate_bins as _ab
    _src = inspect.getsource(_ab.generate_annotated_fasta)
    _old = ('            elif annotation["rank"] == "D":\n'
            '                annotation_str += "; %s (db=%s)" % (annotation.pfam_hits, "pfam")')
    _new = ('            elif annotation["rank"] == "D":\n'
            '                if "pfam_hits" in annotation:\n'
            '                    annotation_str += "; %s (db=%s)" % (annotation.pfam_hits, "pfam")')
    if _old in _src:
        exec(compile(_src.replace(_old, _new), _ab.__file__, "exec"), _ab.__dict__)
except Exception:
    pass
PYEOF

export PYTHONPATH="$PWD:${PYTHONPATH:-}"
exec python3 ./.dram_local.py "$@"
