#!/usr/bin/env python3
"""
Build a valid DRAM CONFIG for annotate/distill from a staged DRAM database
directory, writing it to a local JSON file the caller points DRAM_CONFIG_LOCATION at.

Why this exists: the databases in this project were set up on another host and the
resulting DRAM CONFIG does not travel with the database directory (it lives inside
the setup host's mag_annotator package). Inside the container DRAM therefore falls
back to its packaged default CONFIG, in which every database path is unset — so
annotate silently runs against no databases and distill has no ko_id to summarise.

This locates each database file in the staged directory (by prefix/suffix, so it
is robust to the version-dated filenames and to DRAM's own misspelling of the ETC
sheet as 'etc_mdoule_database') and writes them into a config with a recognised
dram_version so DRAM's load_config() takes the modern path and honours it verbatim.

uniref90 is intentionally left unset (see the pipeline decision): it is by far the
slowest search and contributes only gene descriptions, not KEGG-module coverage,
which is what the distillate needs. viral/vogdb/kegg are likewise left unset.

Usage: build_dram_config.py <dram_db_dir> [output_config.json]
"""
import json
import os
import sys

# (config target, [ (prefix, suffix), ... ]) — first matching file in the DB dir
# wins. The suffixes are chosen precisely so they select the base database file
# and never its mmseqs/hmmer auxiliary siblings, which end in different suffixes
# (e.g. pfam.mmspro vs pfam.mmspro.index/_h/.idx; kofam_profiles.hmm vs .hmm.h3f;
# peptidases.<date>.mmsdb vs .mmsdb.index/.mmsdb_h/.mmsdb.idx).
SEARCH_DATABASES = {
    'kofam_hmm':      [('kofam_profiles', '.hmm')],
    'kofam_ko_list':  [('kofam_ko_list', '.tsv')],
    'dbcan':          [('dbCAN-HMMdb', '.txt')],
    'peptidase':      [('peptidases', '.mmsdb')],
}
DATABASE_DESCRIPTIONS = {
    'pfam_hmm':             [('Pfam-A', '.dat.gz')],
    'dbcan_fam_activities': [('CAZyDB', 'fam-activities.txt')],
    'dbcan_subfam_ec':      [('CAZyDB', 'subfam.ec.txt')],
}
DRAM_SHEETS = {
    'genome_summary_form':   [('genome_summary_form', '.tsv')],
    'module_step_form':      [('module_step_form', '.tsv')],
    'function_heatmap_form': [('function_heatmap_form', '.tsv')],
    'amg_database':          [('amg_database', '.tsv')],
    # DRAM's own setup writes this misspelled; the code reads the correct key.
    'etc_module_database':   [('etc_module_database', '.tsv'),
                              ('etc_mdoule_database', '.tsv')],
}
# Databases deliberately left unset:
#  - uniref:  by far the slowest search; contributes only gene descriptions,
#             not KEGG-module coverage, so not needed for the distillate.
#  - pfam:    the mmseqs *profile* search on this DB's pfam.mmspro crashes the
#             mmseqs build in the DRAM 1.4.6 container ("Score of forward/backward
#             SW differ"). It is not needed for the distillate (module coverage
#             comes from kofam), so it is excluded rather than block the run.
#  - viral/vogdb/kegg: not used for MAG metabolism annotation here.
UNSET_SEARCH = ['uniref', 'pfam', 'viral', 'vogdb', 'kegg']


def pick(dbdir, prefixes):
    """Return abspath of the first file in dbdir matching any (prefix, suffix)."""
    try:
        names = sorted(os.listdir(dbdir))
    except OSError:
        return None
    for prefix, suffix in prefixes:
        for fn in names:
            if fn.startswith(prefix) and fn.endswith(suffix):
                return os.path.abspath(os.path.join(dbdir, fn))
    return None


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("usage: build_dram_config.py <dram_db_dir> [out.json]\n")
        return 2
    dbdir = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else 'LOCAL_DRAM_CONFIG.json'

    # Seed structure from the packaged default so every expected key exists.
    cfg = {}
    try:
        from pkg_resources import resource_filename
        default = resource_filename('mag_annotator', 'CONFIG')
        if os.path.isfile(default):
            with open(default) as fh:
                cfg = json.load(fh)
    except Exception as e:
        sys.stderr.write("WARN: could not load packaged default CONFIG: %s\n" % e)
    if not isinstance(cfg, dict):
        cfg = {}

    sd = cfg.setdefault('search_databases', {})
    dd = cfg.setdefault('database_descriptions', {})
    sh = cfg.setdefault('dram_sheets', {})

    for key, prefixes in SEARCH_DATABASES.items():
        sd[key] = pick(dbdir, prefixes)
    for key, prefixes in DATABASE_DESCRIPTIONS.items():
        dd[key] = pick(dbdir, prefixes)
    for key, prefixes in DRAM_SHEETS.items():
        sh[key] = pick(dbdir, prefixes)
    for key in UNSET_SEARCH:
        sd[key] = None

    desc = pick(dbdir, [('description_db', '.sqlite')])
    cfg['description_db'] = desc

    # A recognised dram_version keeps load_config() on the modern path, which
    # preserves the dicts above instead of treating the file as pre-1.4.0.
    if not cfg.get('dram_version'):
        cfg['dram_version'] = '1.4.0'

    with open(out, 'w') as fh:
        json.dump(cfg, fh, indent=2)

    sys.stderr.write("DRAM config written to %s:\n" % out)
    for label, d, keys in (('search_databases', sd, list(SEARCH_DATABASES) + UNSET_SEARCH),
                           ('database_descriptions', dd, list(DATABASE_DESCRIPTIONS)),
                           ('dram_sheets', sh, list(DRAM_SHEETS))):
        for k in keys:
            sys.stderr.write("  %s.%s -> %s\n" % (label, k, d.get(k)))
    sys.stderr.write("  description_db -> %s\n" % cfg['description_db'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
