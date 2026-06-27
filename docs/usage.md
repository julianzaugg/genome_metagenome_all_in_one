# Usage

## Invocation

```bash
nextflow run . -profile <profile> --mode <mode> --input samplesheet.csv --outdir results
```

- `--mode` — one of `illumina_metagenome`, `nanopore_metagenome`,
  `illumina_isolate`, `nanopore_isolate`, or `download_dbs`.
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
| `host_ref`   | optional                                | leave blank for now; Illumina metagenome uses global `--host_ref` / `--cleanifier_db` |

Per-mode required-input rules are enforced at runtime in `subworkflows/local/input_check.nf`.

For isolate workflows, the samplesheet is the sequencing inventory. Comparative
membership can be controlled separately with `--comparison_manifest`; when that
file is omitted, samples are compared by samplesheet `group` and no external
references are included.

### Hybrid (long + short) reads

A Nanopore-mode sample may **also** provide `fastq_1`/`fastq_2`. Short-read-dependent
steps (e.g. polypolish polishing) then run for that sample only; samples without
short reads skip them. A mixed samplesheet works in a single run — gating is
per-sample via `meta.has_short_reads`.

Dorado basecalling is optional and row-scoped. The workflow uses `long_reads`
FASTQ when present. It schedules Dorado basecalling only for rows with `pod5_dir`
and no `long_reads`, or for all rows with `pod5_dir` when
`--force_dorado_basecalling true`. Dorado polishing is controlled separately by
`--skip_dorado_polish`.

Example (`assets/samplesheets/nanopore_isolate.csv`):

```csv
sample,group,fastq_1,fastq_2,long_reads,pod5_dir,host_ref
ISO_1,speciesX,,,reads/iso1_long.fastq.gz,,
ISO_2,speciesX,reads/iso2_R1.fastq.gz,reads/iso2_R2.fastq.gz,reads/iso2_long.fastq.gz,,
```

### Isolate comparison manifest

Use `--comparison_manifest manifest.csv` when you need mixed species,
overlapping groups, or external reference genomes. Columns:

| Column | Required | Notes |
|--------|----------|-------|
| `comparison_group` | yes | group name; the same sample/reference can appear in multiple groups |
| `entry_type` | yes | `sample` or `reference` |
| `id` | yes | samplesheet sample id, or reference id |
| `fasta` | reference rows | reference genome FASTA |
| `gff` | optional | reference GFF3 for Panaroo |
| `faa` | optional | reference proteins for AMRFinderPlus/chewBBACA context |
| `parsnp_reference` | optional | at most one true row per group |

If no `parsnp_reference=true` row is present for a group, the workflow chooses
the first reference, then the first sample. Reference annotation input is
de-duplicated by `id`, then expanded into each requested group. chewBBACA input
genome names longer than its practical limit are hashed, and the mapping is
published with the chewBBACA outputs.

## Turning steps on/off

Steps are toggled with `--skip_*` / `--run_*` params, e.g.
`--skip_singlem --skip_mobile_elements --run_nonpareil false`. See
`nextflow run . --help` for the full list. Reorderable steps are independent
subworkflows fed from shared upstream channels, so changing the flow is a wiring
edit in `workflows/<mode>.nf`, not a rewrite.

Host removal is on by default for metagenomes. Keep it on for normal runs and
provide either `--cleanifier_db` or `--host_ref`; use `--skip_host_removal true`
only when you deliberately want assembly from QC'd, unfiltered reads.

RPKM is on by default for Illumina metagenomes and requires fastp, assembly, and
the gene catalogue. It uses one selected R1 stream: fastp + host-filtered reads
when host removal runs, otherwise fastp reads. Tune the DIAMOND input length
filter with `--rpkm_min_read_length` or disable the stage with `--skip_rpkm true`.

Raw gzipped FASTQs are checked before QC/read profiling with `gzip -t` and
`seqkit stats`. SingleM `Unexpected line format for DIAMOND output line` errors
can indicate corrupt or malformed FASTQ input, so fix the source reads and resume.

## Tool-choice slots

`--catalogue_clusterer {cdhit|mmseqs}`, `--core_alignment {parsnp|snippy}`,
`--tree_builder {iqtree|raxml}`.

## Resuming

Add `-resume`. Note: Aviary and autocycler are single coarse-grained processes
(they wrap their own pipelines), so the whole process is the resume unit — their
internal steps are opaque to Nextflow.
