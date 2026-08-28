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

### Comparison samplesheets (`illumina_metagenome` only)

Two optional, independent samplesheets let you compare a second dataset against
this run's own results without assembling or binning it — see
[output.md](output.md#comparison-reads-and-comparison-assemblies-29_comparison_reads-30_comparison_assemblies)
for what each produces.

`--comparison_reads reads.csv` — one row per external sample to QC/host-remove and
map:

| Column     | Required | Notes |
|------------|----------|-------|
| `sample`   | yes      | unique id — must not collide with `--input` or `--comparison_assemblies` |
| `fastq_1`  | yes      | gzipped FASTQ |
| `fastq_2`  | yes      | gzipped FASTQ |

`--comparison_assemblies assemblies.csv` — one row per external, pre-binned assembly
to predict genes from:

| Column     | Required | Notes |
|------------|----------|-------|
| `sample`   | yes      | unique id — must not collide with `--input` or `--comparison_reads` |
| `assembly` | yes      | assembly FASTA (contigs/scaffolds), optionally gzipped |

Either can be used alone, and rows in one need not correspond to rows in the other
or in the main samplesheet.

## Turning steps on/off

Steps are toggled with `--skip_*` / `--run_*` params, e.g.
`--skip_singlem --skip_mobile_elements --run_nonpareil false`. See
[parameters.md](parameters.md) for the full list, broken down by which
`--mode` each flag applies to and the requirements/conflicts between them.
Reorderable steps are independent subworkflows fed from shared upstream
channels, so changing the flow is a wiring edit in `workflows/<mode>.nf`, not
a rewrite.

Host removal is on by default for metagenomes. Keep it on for normal runs and
provide either `--cleanifier_db` or `--host_ref`; use `--skip_host_removal true`
only when you deliberately want assembly from QC'd, unfiltered reads.

RPKM is on by default for Illumina metagenomes and requires fastp, assembly, and
the gene catalogue. It uses one selected R1 stream: fastp + host-filtered reads
when host removal runs, otherwise fastp reads. Tune the DIAMOND input length
filter with `--rpkm_min_read_length` or disable the stage with `--skip_rpkm true`.
If available, provide the prebuilt SingleM marker DIAMOND directory and marker
stats TSV with `--rpkm_singlem_marker_dbs` and
`--rpkm_singlem_marker_lengths`.

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

### Pin the session on long runs

A bare `-resume` resumes **the most recent session in that directory**, which is
not necessarily your last real run. Every `nextflow run` opens a session,
including ones that execute nothing:

- `nextflow run . --help` — `--help` is a pipeline *parameter*, so this is still a run
- any launch that dies during config/schema validation

Either will become "the most recent session", and the next bare `-resume` then
attaches to it, finds an empty cache, and silently re-runs the pipeline from the
first process. The symptom — early steps like `FASTQ_GZIP_TEST` and `FASTP`
re-running right after a `git pull` — looks like a code regression but is not.

Check before resuming anything expensive, and pass the session explicitly:

```bash
nextflow log | tail -5          # confirm the SESSION ID of your last real run
nextflow run . ... -resume <session-id>
```

To inspect parameters without opening a session, read [parameters.md](parameters.md)
or run `nextflow config .` — neither is a `nextflow run`.

`DUMP_SOFTWARE_VERSIONS` re-runs whenever modules are added or removed. That one
is expected and does not indicate a broken cache.
