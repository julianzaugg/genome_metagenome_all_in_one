# Reference databases

Each database is a parameter, defaulting to the `/srv/db` layout already in use.
Override per profile (Bunya paths differ from `/srv`) or on the command line.

| Param                   | Tool(s)            | Default (`/srv`) |
|-------------------------|--------------------|------------------|
| `--gtdbtk_db`           | GTDB-Tk, Aviary    | `/srv/db/gtdbtk/official/release232` |
| `--checkm2_db`          | CheckM2, Aviary    | `/srv/db/checkm2_data/1.0.2/CheckM2_database` |
| `--checkm1_db`          | CheckM1            | `/srv/db/checkm1` |
| `--singlem_metapackage` | SingleM            | `/srv/db/singlem_packages/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb` |
| `--sylph_db`            | sylph              | path/glob of the `.syldb` file(s) to use (see note) |
| `--bakta_db`            | Bakta              | `/srv/db/bakta/6.0/db` |
| `--dram_db`             | DRAM               | `/srv/db/DRAM_1.4.6` |
| `--genomad_db`          | geNomad            | `/srv/db/genomad` |
| `--checkv_db`           | CheckV             | `/srv/db/checkv/checkv-db-v1.5` |
| `--eggnog_db`           | Aviary             | `/srv/db/eggnog/5.0` |
| `--amrfinder_db`        | AMRFinderPlus      | `null` → container's bundled DB |
| `--host_ref`            | Cleanifier         | host FASTA used to build an index if `--cleanifier_db` is unset |
| `--cleanifier_db`       | Cleanifier         | `null` (use a prebuilt `.filter`, otherwise build from `--host_ref`) |
| `--cleanifier_nobjects` | Cleanifier         | optional `cleanifier index -n` override; defaults to estimated FASTA bases |
| `--genomespot_models`   | GenomeSPOT         | `null` |

### Cleanifier host removal

Host removal is enabled by default. Provide one of:

```bash
--cleanifier_db /path/to/prebuilt.cleanifier.filter
# or
--host_ref /path/to/host.fa.gz
```

When `--host_ref` is used, the pipeline first builds a Cleanifier index and
publishes it under `04_host_removed/index/`. For repeated runs, point
`--cleanifier_db` at the published `.filter` file to avoid rebuilding. Cleanifier
requires `-n/--nobjects` while indexing; by default the pipeline estimates this
from FASTA base count. Indexing uses `-k 29 --windowsize 33` by default; adjust
`CLEANIFIER_INDEX.ext.args` in `conf/modules.config` if you need a different
mask/k-mer setup. For `.zst` FASTA or tighter sizing, set `--cleanifier_nobjects`
explicitly.

### sylph `.syldb` selection

`--sylph_db` is a path **or glob** of the `.syldb` file(s) to profile against —
you choose which, because a sylph DB folder typically holds several (e.g.
`/srv/db/sylph` has `gtdb-r226`, `gtdb-r232` and `fungi`). The pipeline profiles
against **exactly the files the glob matches** (it does not blindly use the whole
folder). Examples:

```bash
--sylph_db '/srv/db/sylph/gtdb-r232-c200-dbv1.syldb'        # one DB
--sylph_db '/srv/db/sylph/{gtdb-r232,fungi}*.syldb'         # GTDB r232 + fungi (not r226)
```

### sylph-tax taxonomy

Set **`--sylph_tax_metadata`** to metadata matching the `.syldb` files selected
by `--sylph_db` to run the taxonomy chain after profiling (`sylph-tax taxprof`
→ `sylph-tax merge`), producing `merged_relative_abundance.tsv` and
`merged_sequence_abundance.tsv` in `02_sylph/`. Leave it `null` to run
`sylph profile` only.

```bash
--sylph_tax_metadata '/srv/db/sylph/gtdb_r232_metadata.tsv.gz'
```
For multi-DB runs, use a glob or list that covers exactly the selected databases,
for example GTDB r232 plus fungi. `conf/local.config` pre-sets this to
`/srv/db/sylph/gtdb_r232_metadata.tsv.gz`.

Tools that read DB locations from environment variables
(`GTDBTK_DATA_PATH`, `CHECKM2_DATA_PATH`, `SINGLEM_METAPACKAGE_PATH`,
`EGGNOG_DATA_DIR`) have them exported inside the module from the corresponding
param, so the param is the single source of truth.

## Downloading databases

You can point the params at existing locations (no download), **or** build them:

```bash
nextflow run . -profile bunya --mode download_dbs --download_db all --outdir refs
# or a subset:
nextflow run . -profile bunya --mode download_dbs --download_db gtdbtk,checkv,genomad
```

Databases are written under `<outdir>/<db_outdir>/<name>` (`--db_outdir` defaults
to `databases`). Then set e.g. `--gtdbtk_db refs/databases/gtdbtk`.

The downloader commands in `modules/local/download_dbs.nf` are best-effort for
the pinned tool versions — confirm the exact flags for your versions.
