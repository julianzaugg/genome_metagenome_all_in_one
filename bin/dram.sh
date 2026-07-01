#!/usr/bin/env bash
#
# Wrapper around DRAM.py that works around a DRAM 1.5.0 CLI bug.
#
# DRAM 1.5.0's annotate / annotate_genes / distill subparsers define a
# --config_loc argument (default None), but its __main__ calls the handler with
# every parsed arg except 'func':
#     args_dict = {i: j for i, j in vars(args).items() if i != 'func'}
#     args.func(**args_dict)
# The handler functions (e.g. annotate_called_genes_cmd) don't accept config_loc,
# so this raises `TypeError: ... unexpected keyword argument 'config_loc'` on
# every invocation — regardless of whether --config_loc is passed.
#
# We run a local copy of DRAM.py with config_loc also stripped from that dict.
# We supply the config via the DRAM_CONFIG_LOCATION env var instead, so dropping
# the flag loses nothing. On DRAM 1.4.6 the sed simply doesn't match and the
# copy runs unchanged, so this wrapper is safe across both versions.
set -euo pipefail

src=$(command -v DRAM.py)
cp "$src" ./.dram_local.py
sed -i "s/if i != 'func'/if i not in ('func', 'config_loc')/" ./.dram_local.py
exec python3 ./.dram_local.py "$@"
