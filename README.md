# gmaio — genome & metagenome all-in-one

A single Nextflow (DSL2) pipeline consolidating four analysis tracks:

| `--mode`               | Platform | Sample type | Status (v1)        |
|------------------------|----------|-------------|--------------------|
| `illumina_metagenome`  | Illumina | metagenome  | **built end-to-end** |
| `nanopore_metagenome`  | Nanopore | metagenome  | scaffold (stubs)   |
| `illumina_isolate`     | Illumina | isolate     | scaffold (stubs)   |
| `nanopore_isolate`     | Nanopore | isolate     | scaffold (stubs)   |

It is **nf-core-style but standalone** — modules/subworkflows layout, samplesheet
schema, `nf-schema` param validation — without nf-core governance. Reuses
nf-core/modules where they exist; bespoke tools (Aviary, autocycler, myloasm,
dorado, sylph, cleanifier, DRAM, chewbacca, …) are custom local modules.

**Containers only.** Reused nf-core modules pull biocontainers from a registry;
bespoke tools use local `.sif` images on a shared path. See
[docs/containers.md](docs/containers.md).

## Quick start

```bash
# Stub dry-run (no real tools, validates wiring)
nextflow run . -profile test,local --mode illumina_metagenome \
  --input assets/samplesheets/illumina_metagenome.csv -stub

# Real run on Bunya
nextflow run . -profile bunya --mode illumina_metagenome \
  --input samplesheet.csv --outdir results
```

## Profiles

- `bunya` — Bunya HPC, SLURM executor, account `a_ace` (primary target).
- `bunya_gpu` — adds H100 GPU resources for dorado steps (`gpu_cuda`, `gpu` qos).
- `local` — the `/srv` Linux server (local executor).
- `test` — tiny resources + bundled mini data for development.

## Documentation

- [docs/configuration.md](docs/configuration.md) — what each config file is and what to edit for setup.
- [docs/usage.md](docs/usage.md) — per-mode usage and the samplesheet spec (incl. hybrid long+short).
- [docs/databases.md](docs/databases.md) — DB params and the optional download workflow.
- [docs/containers.md](docs/containers.md) — `.sif` image checklist.
- [docs/output.md](docs/output.md) — output layout (mirrors the numbered-dir convention).

See [CLAUDE.md](CLAUDE.md) for a concise repo map and conventions.
