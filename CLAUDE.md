# CLAUDE.md

Concise guide for working in this repo. Keep this file minimal.

## What this is
`gmaio` — one Nextflow DSL2 pipeline, four `--mode` tracks:
`illumina_metagenome` (built), `nanopore_metagenome`, `illumina_isolate`,
`nanopore_isolate` (scaffolds). nf-core-*style* but standalone. Containers only.

## Repo map
- `main.nf` — entry; routes on `params.mode` to a `workflows/*.nf`.
- `workflows/` — one file per mode.
- `subworkflows/local/` — reusable step groups (qc, profiling, assembly, binning, …).
- `modules/local/` — bespoke tool wrappers (Aviary, autocycler, myloasm, dorado, sylph, cleanifier, DRAM, chewbacca, checkv_cluster, …).
- `modules/nf-core/` — installed via `nf-core modules install`.
- `bin/` — helper scripts (python/R), on `PATH` inside processes.
- `conf/` — `base`, `modules` (ext.args + publishDir), `containers`, `bunya`, `bunya_gpu`, `local`, `test`.
- `assets/schema_input.json` — samplesheet validation. `nextflow_schema.json` — params.

## Run
Use Nextflow **25.04.x** (matches Bunya and the Linux server). Pin locally with
`export NXF_VER=25.04.6`. The config is already migrated to the modern idioms
(`process.resourceLimits`, no `new Date()`), so it parses under Nextflow **26+**'s
strict config parser. The **scripts** are not yet migrated, though — 26+'s strict
script parser rejects e.g. the `switch` in `main.nf` — so on 26+ run with
`NXF_SYNTAX_PARSER=v1` until the `.nf` files are migrated too.
```bash
# wiring check, no real tools (needs DB params set, e.g. on a server)
nextflow run . -profile test,local --mode illumina_metagenome \
  --input assets/samplesheets/illumina_metagenome.csv -stub
# laptop check without databases: skip the steps that require reference DBs
nextflow run . -profile test,local --mode illumina_metagenome \
  --input assets/samplesheets/illumina_metagenome.csv -stub \
  --skip_host_removal true --skip_assembly true --skip_binning true \
  --skip_gene_catalogue true --skip_mobile_elements true --skip_rpkm true
# real
nextflow run . -profile bunya --mode <mode> --input samplesheet.csv --outdir results
```

## Output directory numbering (`conf/modules.config`)
All `publishDir` paths use a `NN_name` prefix. Numbers must be unique **within
each mode** (metagenome and isolate run separately, so they can reuse numbers).
After adding or renaming any `publishDir`, verify no prefix collision exists:
```bash
grep -oP '\d{2,3}_\w+' conf/modules.config | sort -V | uniq -c | sort -rn | awk '$1>1'
```
A non-empty result means a collision — renumber to resolve before committing.
Update `docs/output.md` to match whenever numbers change.

## Conventions
- Tool flags live in `conf/modules.config` as `ext.args`, never hard-coded in process bodies.
- Containers set per-process in `conf/containers.config` (biocontainer URI or local `.sif`).
- Steps are toggled with `--skip_*` / `--run_*` params; subworkflows guarded by `if`.
- Cross-sample steps (dereplication, GTDB-Tk, gene catalogue, pangenome, ANI) use
  `.collect()`; group-scoped comparisons use `groupTuple()` on the `group` column.
- Hybrid: a Nanopore sample may also have `fastq_1/2`; `meta.has_short_reads` gates
  short-read-dependent steps (e.g. polypolish) per-sample.

## Setup before a real run (see docs/)
- Provide `.sif` paths for bespoke tools — `docs/containers.md`.
- Set DB params (GTDB-Tk, CheckM2, singlem, bakta, DRAM, geNomad, CheckV, host) or
  run `--mode download_dbs` — `docs/databases.md`.
- Confirm Bunya account/partition/qos in `conf/bunya*.config`.
