# Usage

## Invocation

```bash
nextflow run . -profile <profile> --mode <mode> --input samplesheet.csv --outdir results
```

- `--mode` — one of `illumina_metagenome` (built), `nanopore_metagenome`,
  `illumina_isolate`, `nanopore_isolate` (scaffolds), or `download_dbs`.
- `-profile` — `bunya` (primary), `bunya_gpu` (adds H100 for dorado), `local`
  (the `/srv` server), or `test` (tiny, self-contained, for `-stub`).
- `--input` — samplesheet (below). Not needed for `download_dbs`.

Run `nextflow run . --help` to see all parameters (rendered from `nextflow_schema.json`).

## Samplesheet

One row per sample. Columns:

| Column       | Required when…                          | Notes |
|--------------|-----------------------------------------|-------|
| `sample`     | always                                  | unique sample id |
| `group`      | optional                                | keys cross-sample comparison sets (pangenome / parsnp / ANI). Defaults to `all`. |
| `fastq_1`    | Illumina modes; optional for Nanopore   | gzipped FASTQ |
| `fastq_2`    | with `fastq_1`                          | gzipped FASTQ |
| `long_reads` | Nanopore modes (or `pod5_dir`)          | gzipped FASTQ |
| `pod5_dir`   | Nanopore modes (alternative to `long_reads`) | directory of POD5 for dorado basecalling |
| `host_ref`   | optional                                | per-sample host reference FASTA |

Per-mode required-input rules are enforced at runtime in `subworkflows/local/input_check.nf`.

### Hybrid (long + short) reads

A Nanopore-mode sample may **also** provide `fastq_1`/`fastq_2`. Short-read-dependent
steps (e.g. polypolish polishing) then run for that sample only; samples without
short reads skip them. A mixed samplesheet works in a single run — gating is
per-sample via `meta.has_short_reads`.

Example (`assets/samplesheets/nanopore_isolate.csv`):

```csv
sample,group,fastq_1,fastq_2,long_reads,pod5_dir,host_ref
ISO_1,speciesX,,,reads/iso1_long.fastq.gz,,
ISO_2,speciesX,reads/iso2_R1.fastq.gz,reads/iso2_R2.fastq.gz,reads/iso2_long.fastq.gz,,
```

## Turning steps on/off

Steps are toggled with `--skip_*` / `--run_*` params, e.g.
`--skip_singlem --skip_mobile_elements --run_nonpareil false`. See
`nextflow run . --help` for the full list. Reorderable steps are independent
subworkflows fed from shared upstream channels, so changing the flow is a wiring
edit in `workflows/<mode>.nf`, not a rewrite.

## Tool-choice slots

`--catalogue_clusterer {cdhit|mmseqs}`, `--core_alignment {parsnp|snippy}`,
`--tree_builder {iqtree|raxml}`.

## Resuming

Add `-resume`. Note: Aviary and autocycler are single coarse-grained processes
(they wrap their own pipelines), so the whole process is the resume unit — their
internal steps are opaque to Nextflow.
