# Configuration — what each file is, and what to edit

How settings combine: **`nextflow.config` defaults → the profile you select
(`conf/<profile>.config`) → `--flags` on the command line**. Later wins. So you
rarely edit the defaults — you put machine-specific settings in a profile and
one-off settings on the command line.

## The files

| File | What it controls | Do you edit it? |
|------|------------------|-----------------|
| **`nextflow.config`** | Top level: declares every parameter + its **default**, loads the `conf/*` files, wires the profiles, turns Apptainer on. | Rarely — only to change a *global default* or add a new param. |
| **`nextflow_schema.json`** | Validation + `--help` text for parameters. | Only when you **add/rename a param** (keep it in sync with `nextflow.config`). |
| **`conf/base.config`** | Resource **labels** — how many CPUs/memory/time each `process_low/medium/high/maxmem/gpu` tier gets, and retry behaviour. | To tune resources for your hardware. |
| **`conf/modules.config`** | Per-process **tool flags** (`ext.args`) and **output locations** (`publishDir`, the numbered dirs). | To change a tool's flags or where outputs land. **This is where tool options live — never hard-code them in a module.** |
| **`conf/containers.config`** | Per-process **container image** (quay.io URI or local `.sif`). | To change an image version/source, or point a bespoke tool at your `.sif`. |
| **`conf/local.config`** | The **`local` profile** — your `/srv` server: executor, CPU/mem ceilings, project-local container/cache paths, and **database paths**. | **Yes — this is your main setup file for the local server.** |
| **`conf/bunya.config`** / **`conf/bunya_gpu.config`** | The **`bunya`** profiles — SLURM account/partition/qos, project-local container/cache paths, GPU options. | Yes, if you run on Bunya (confirm account/qos). |
| **`conf/test.config`** | The **`test`** profile — tiny resources + stub data for `-stub` wiring checks. | No (leave as-is). |
| **`assets/schema_input.json`** | Allowed **samplesheet** columns + validation. | Only if you change the samplesheet format. |
| **`docs/*.md`** | Documentation only. | Editing these changes **nothing** in the pipeline. |

## "I want to change X" → edit this

| I want to… | Edit |
|------------|------|
| Set where my **databases** are | profile config (`conf/local.config` or a Bunya overlay) — or `--<db>_db` on the CLI |
| Point a tool at **my own `.sif`** | `conf/containers.config` (the `withName` line) |
| Use a **different tool version** (container) | `conf/containers.config` (change the tag/URI) |
| Change a **tool's flags** (e.g. cd-hit identity, fastp settings) | `conf/modules.config` (`ext.args`) |
| Change an **output folder name/number** | `conf/modules.config` (`publishDir`) |
| Give a step **more CPU/RAM/time** | `conf/base.config` (its `withLabel` tier) — or `--max_cpus/--max_memory` |
| Turn a **step on/off** | a `--skip_*` / `--run_*` flag (see `nextflow run . --help`) |
| Set **SLURM account/partition** | `conf/bunya.config` |
| Where bespoke **`.sif` images live** | `conf/local.config` / `conf/bunya.config` (`container_base`, default `<projectDir>/containers`) |

## Minimum setup to run on the local server

1. **`conf/local.config`** — confirm the DB paths (pre-filled to `/srv/db/...`)
   and the CPU/mem ceilings. Container images default to project-local
   `containers/` and `.apptainer_cache/`.
2. **Host removal** — because `--skip_host_removal false` is the default, set
   either `cleanifier_db` to a prebuilt `.filter` index or `host_ref` to a FASTA
   that Cleanifier can index. For one-off runs, pass `--host_ref` or
   `--cleanifier_db` on the command line instead.
3. **Containers** — most tools auto-pull from quay.io (nothing to do). Build only
   the local `.sif` images if you need them (`dorado`, `genomespot`) — see
   [containers.md](containers.md).
4. Run:

```bash
nextflow run . -profile local --mode illumina_metagenome \
  --input s.csv --outdir results \
  --host_ref /srv/db/host/host.fa.gz
```

If you already built a Cleanifier index, use `--cleanifier_db /path/to/index.filter`
instead of `--host_ref`.

## Minimum setup to run on Bunya

1. **`conf/bunya.config`** — confirm `--account`, queue/partition and resource
   ceilings for your allocation. Container images default to project-local
   `containers/` and `.apptainer_cache/`.
2. **Database paths** — Bunya jobs cannot use `/srv/db` unless that path is mounted
   on the cluster. Pass Bunya-visible paths on the command line, or create a small
   site/profile config that overrides the DB params from [databases.md](databases.md).
3. **Host removal** — pass either a Bunya-visible `--cleanifier_db` prebuilt index
   or a Bunya-visible `--host_ref` FASTA. Building the index from FASTA can be
   expensive; reuse a published `.filter` for repeated runs.
4. **Work/output locations** — put `--outdir` on project/scratch storage, not a
   login-node home directory.
5. Run:

```bash
nextflow run . -profile bunya --mode illumina_metagenome \
  --input s.csv \
  --outdir /scratch/project/a_ace/$USER/gmaio_results \
  --host_ref /scratch/project/a_ace/db/host/host.fa.gz \
  --gtdbtk_db /scratch/project/a_ace/db/gtdbtk/official/release226 \
  --checkm2_db /scratch/project/a_ace/db/checkm2_data/1.0.2/CheckM2_database \
  --singlem_metapackage /scratch/project/a_ace/db/singlem_packages/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb \
  --sylph_db /scratch/project/a_ace/db/sylph/gtdb-r232-c200-dbv1.syldb \
  --dram_db /scratch/project/a_ace/db/DRAM_1.4.6
```

Add the other DB overrides from [databases.md](databases.md) if their defaults are
not visible on Bunya.

Anything you want different for a single run, append as `--flag value` — no file
edit needed.
