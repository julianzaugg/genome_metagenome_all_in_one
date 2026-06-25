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
| **`conf/local.config`** | The **`local` profile** — your `/srv` server: executor, CPU/mem ceilings, `container_base`, Apptainer cache, and **database paths**. | **Yes — this is your main setup file for the local server.** |
| **`conf/bunya.config`** / **`conf/bunya_gpu.config`** | The **`bunya`** profiles — SLURM account/partition/qos, cache dir, GPU options. | Yes, if you run on Bunya (confirm account/qos). |
| **`conf/test.config`** | The **`test`** profile — tiny resources + stub data for `-stub` wiring checks. | No (leave as-is). |
| **`assets/schema_input.json`** | Allowed **samplesheet** columns + validation. | Only if you change the samplesheet format. |
| **`docs/*.md`** | Documentation only. | Editing these changes **nothing** in the pipeline. |

## "I want to change X" → edit this

| I want to… | Edit |
|------------|------|
| Set where my **databases** are | `conf/local.config` (`params { … }`) — or `--<db>_db` on the CLI |
| Point a tool at **my own `.sif`** | `conf/containers.config` (the `withName` line) |
| Use a **different tool version** (container) | `conf/containers.config` (change the tag/URI) |
| Change a **tool's flags** (e.g. cd-hit identity, fastp settings) | `conf/modules.config` (`ext.args`) |
| Change an **output folder name/number** | `conf/modules.config` (`publishDir`) |
| Give a step **more CPU/RAM/time** | `conf/base.config` (its `withLabel` tier) — or `--max_cpus/--max_memory` |
| Turn a **step on/off** | a `--skip_*` / `--run_*` flag (see `nextflow run . --help`) |
| Set **SLURM account/partition** | `conf/bunya.config` |
| Where bespoke **`.sif` images live** | `conf/local.config` / `conf/bunya.config` (`container_base`) |

## Minimum setup to run on the local server

1. **`conf/local.config`** — confirm the DB paths (pre-filled to `/srv/db/...`),
   `container_base`, `apptainer.cacheDir`, and the CPU/mem ceilings.
2. **Containers** — most tools auto-pull from quay.io (nothing to do). Build only
   the 3 local `.sif` images if you need them (`dorado`, `genomespot`,
   `minimap2_samtools`) — see [containers.md](containers.md).
3. Run: `nextflow run . -profile local --mode illumina_metagenome --input s.csv --outdir results`.

Anything you want different for a single run, append as `--flag value` — no file
edit needed.
