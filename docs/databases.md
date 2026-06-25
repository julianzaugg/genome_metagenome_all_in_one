# Reference databases

Each database is a parameter, defaulting to the `/srv/db` layout already in use.
Override per profile (Bunya paths differ from `/srv`) or on the command line.

| Param                   | Tool(s)            | Default (`/srv`) |
|-------------------------|--------------------|------------------|
| `--gtdbtk_db`           | GTDB-Tk, Aviary    | `/srv/db/gtdbtk/official/release226` |
| `--checkm2_db`          | CheckM2, Aviary    | `/srv/db/checkm2_data/1.0.2/CheckM2_database` |
| `--checkm1_db`          | CheckM1            | `/srv/db/checkm1` |
| `--singlem_metapackage` | SingleM            | `/srv/db/singlem_packages/S5.4.0...smpkg.zb` |
| `--sylph_db`            | sylph              | `/srv/db/sylph/gtdb-r226` (dir of `*.syldb`) |
| `--bakta_db`            | Bakta              | `/srv/db/bakta/6.0/db` |
| `--dram_db`             | DRAM               | `/srv/db/dram` |
| `--genomad_db`          | geNomad            | `/srv/db/genomad` |
| `--checkv_db`           | CheckV             | `/srv/db/checkv/checkv-db-v1.5` |
| `--eggnog_db`           | Aviary             | `/srv/db/eggnog/5.0` |
| `--amrfinder_db`        | AMRFinderPlus      | `null` → container's bundled DB |
| `--host_ref`            | host removal       | `null` (or per-sample `host_ref` column) |
| `--cleanifier_db`       | cleanifier         | `null` |
| `--genomespot_models`   | GenomeSPOT         | `null` |

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
