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
```bash
# wiring check, no real tools
nextflow run . -profile test,local --mode illumina_metagenome \
  --input assets/samplesheets/illumina_metagenome.csv -stub
# real
nextflow run . -profile bunya --mode <mode> --input samplesheet.csv --outdir results
```

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
- Confirm Bunya account/partition/qos and `apptainer.cacheDir` in `conf/bunya*.config`.
