# Parameters

Every parameter is also documented via `nextflow run . --help` (rendered from
`nextflow_schema.json`). This doc adds the context `--help` can't show: which
of the four `--mode` tracks each flag applies to, what a step toggle
concretely skips or enables, and the requirements/conflicts between flags
that would otherwise surprise you at runtime. For database setup see
[databases.md](databases.md); for invocation/samplesheet basics see
[usage.md](usage.md).

## General / pipeline-level params

| Param | Default | What it does |
|-------|---------|---------------|
| `--mode` | — (required) | Which track to run: `illumina_metagenome`, `nanopore_metagenome`, `illumina_isolate`, `nanopore_isolate`, or `download_dbs`. |
| `--input` | — | Samplesheet CSV path. Not required for `--mode download_dbs`. |
| `--outdir` | `results` | Output directory. |
| `--publish_dir_mode` | `copy` | Nextflow `publishDir` mode (`symlink`, `link`, `copy`, `move`, …). |
| `--container_base` | — | Directory holding local `.sif` images for bespoke tools (Aviary, Dorado, GenomeSPOT). See [containers.md](containers.md). |
| `--aviary_container` | — | Optional explicit Aviary `.sif` path; if unset, `${container_base}/aviary_0.13.0.sif` is used. |
| `--tracs_container` | — | Optional explicit TRACS `.sif` path; if unset the quay.io biocontainer is used. Only needed if `TRACS_*` crashes with exit status 132 (SIGILL) — upstream compiles with `-march=native`, so the published image is tied to bioconda's build-host CPU. See [containers.md](containers.md). |

## Step toggles, by mode

All `--skip_*`/`--run_*` flags are booleans. A `--skip_*` flag defaults to
`false` (the step runs) unless noted; `--run_*` flags list their actual
default.

### Shared across all 4 modes

| Param | Default | What it skips/enables |
|-------|---------|------------------------|
| `--skip_qc` | `false` | FASTP (Illumina) or FASTPLONG (Nanopore, via long-read QC) read trimming. |
| `--skip_assembly` | `false` | Assembly: metaSPAdes (`illumina_metagenome`), Shovill (`illumina_isolate`), myloasm (`nanopore_metagenome`), Autocycler (`nanopore_isolate`). Downstream binning/mapping/gene-catalogue steps get an empty assembly channel if skipped. |
| `--skip_checkm` | `false` | CheckM2 completeness/contamination prediction. Required to be `false` when using `--reference_genomes`. |
| `--skip_taxonomy` | `false` | GTDB-Tk classification (+ name restoration for reference genomes). Must be `false` for `--reference_genomes_taxonomy`, `--marker_tree_include_references`, and `--run_marker_tree` to take effect. |
| `--skip_annotation` | `false` | Metagenome modes: DRAM annotation of the gene catalogue. Isolate modes: Bakta + MLST + AMRFinderPlus + ISEScan together — these four can't be toggled independently in isolate modes. Conflicts with `--skip_comparative false` in isolate modes (comparative analysis needs Bakta's GFF/FAA output). |
| `--skip_mobile_elements` | `false` | The whole mobile-elements subworkflow: geNomad virus/plasmid prediction, CheckV quality assessment, and ANI-based clustering of both. |
| `--skip_read_mapping` | `false` | CoverM contig-level read mapping to sample assemblies (all modes), and in metagenome modes, CoverM genome-level mapping to representative/HQ/dereplicated genome sets. |
| `--skip_mapping_assessment` | `false` | The bases-mapped / percent-of-sequenced-bases metrics (`MAPPING_ASSESS`) computed alongside every CoverM read-mapping step above. Has no effect if `--skip_read_mapping` is set (there is nothing to assess). See `docs/output.md` for the metrics themselves. |
| `--run_checkm1` | `true` | CheckM1 `lineage_wf`, run alongside CheckM2. Feeds `--hq_quality_source` HQ classification. Requires `--checkm1_db`. |
| `--run_genomespot` | `true` | GenomeSPOT genome-trait prediction on representative genomes. |
| `--run_barrnap` | `true` | Barrnap rRNA prediction on representative genomes. |
| `--run_dram_bins` | `false` | Per-bin DRAM annotation + distillation, in addition to gene-catalogue DRAM annotation. |

### Metagenome-only (`illumina_metagenome` + `nanopore_metagenome`)

| Param | Default | What it skips/enables |
|-------|---------|------------------------|
| `--skip_sylph` | `false` | sylph read-level taxonomic profiling of raw reads. |
| `--skip_singlem` | `false` | SingleM community profiling of raw reads. |
| `--skip_host_removal` | `false` | Cleanifier host-read removal. Requires `--cleanifier_db` or `--host_ref` when enabled (errors otherwise). |
| `--skip_binning` | `false` | The entire Aviary binning block — and, because they're nested inside it, everything downstream too: CheckM, dereplication, genome read-mapping, taxonomy/QC, and the marker-gene tree. |
| `--skip_dereplication` | `false` | Cross-sample dereplication (`COVERM_CLUSTER`/`_HQ`/`_HQ_REF`); all bins are used as "representatives" instead. Must be `false` (with `--skip_binning false`) for `--reference_genomes`, `marker_tree_genome_source=hq_representatives`, and `--within_sample_dereplication` other than `none`. |
| `--skip_gene_catalogue` | `false` | Gene prediction (Pyrodigal on scaffolds) and the whole gene-catalogue subworkflow (CD-HIT clustering, CDS extraction, membership tabulation, optional DRAM annotation, and the expanded catalogue when `--reference_genomes` is set). |
| `--skip_rpkm` | `false` | The RPKM subworkflow (SingleM-marker-normalized gene-catalogue abundance via DIAMOND blastx). Requires `--skip_qc false` and either `--skip_assembly false` or `--skip_gene_catalogue false` (hard error otherwise). |
| `--within_sample_dereplication` | `none` | `sample`/`group` enables an independent within-sample or within-group dereplication path, reusing the pooled CheckM reports but clustering only that unit's own bins. Requires `--skip_binning false`. Independent of `--skip_dereplication` — see the schema description for the four across/within combinations. |
| `--run_nonpareil` | `true` | Nonpareil sequencing-coverage/diversity estimation on host-removed clean reads. |
| `--run_marker_tree` | `false` | The marker-gene tree subworkflow (bac120/ar53 marker-protein alignment + closest/related GTDB reference selection + IQTree or VeryFastTree). Requires `--skip_taxonomy false`. See the `marker_tree_*` family below for tuning. |
| `--run_instrain` | `false` | inStrain strain-level comparison of samples against the shared cross-sample dereplicated reference set (bowtie2 -> profile -> compare -> summary), answering whether two samples carry the same strain. Requires `--skip_binning false --skip_dereplication false`. `illumina_metagenome` only. See the `strain_*` family below. |
| `--run_tracs` | `false` | TRACS strain/transmission comparison against the same shared reference set (build-db -> align -> combine -> distance -> cluster), giving pairwise SNP distances and transmission clusters. Builds its reference database from the pipeline's own MAGs, so no GTDB download is needed. Requires `--skip_binning false --skip_dereplication false`. Both metagenome modes; the only option for nanopore. |
| `--reference_genomes_taxonomy` | `false` | Classify `--reference_genomes` with GTDB-Tk alongside the MAGs. Requires `--reference_genomes` and `--skip_taxonomy false`. |
| `--marker_tree_include_references` | `false` | Place `--reference_genomes` in the marker-gene tree alongside the MAGs. Requires `--reference_genomes` and `--skip_taxonomy false`. |

### Isolate-only (`illumina_isolate` + `nanopore_isolate`)

| Param | Default | What it skips/enables |
|-------|---------|------------------------|
| `--skip_comparative` | `false` | The entire isolate-comparative subworkflow (FastANI, ParSNP+Gubbins, Panaroo(+IQTree), chewBBACA) across sample/reference comparison groups. |
| `--run_fastani` | `true` | FastANI within isolate-comparative. Only takes effect when `--skip_comparative false`. |
| `--run_parsnp` | `true` | ParSNP+Gubbins core-genome alignment within isolate-comparative. Only takes effect when `--skip_comparative false`. |
| `--run_panaroo` | `true` | Panaroo pangenome analysis within isolate-comparative. Only takes effect when `--skip_comparative false`. |
| `--run_chewbbaca` | `true` | chewBBACA cgMLST within isolate-comparative. Only takes effect when `--skip_comparative false`. |
| `--run_tree` | `true` | IQTree tree built from the Panaroo core-gene alignment. Only takes effect when `--run_panaroo` is also `true`. |

### Nanopore-only (`nanopore_metagenome` + `nanopore_isolate`)

| Param | Default | What it skips/enables |
|-------|---------|------------------------|
| `--skip_porechop` | `false` | Porechop adapter trimming inside long-read QC. |
| `--skip_dorado_polish` | `false` | Dorado long-read consensus polishing after assembly (Autocycler in isolate mode, myloasm in metagenome mode). |
| `--force_dorado_basecalling` | `false` | Forces POD5 → Dorado basecalling even when direct long-read FASTQ is also supplied for a sample. |

### `nanopore_isolate`-only

| Param | Default | What it skips/enables |
|-------|---------|------------------------|
| `--skip_polypolish` | `false` | Polypolish short-read hybrid polishing of the long-read assembly. Only meaningful for samples that also carry short reads (`meta.has_short_reads`). |

### HQ classification (shared, conditional on CheckM steps)

| Param | Default | What it does |
|-------|---------|---------------|
| `--hq_quality_source` | `both` | Which CheckM report(s) classify a bin as HQ (`checkm1`, `checkm2`, or `both` — a bin is HQ if it passes in any selected report). Only meaningful when `--run_checkm1 true` and/or `--skip_checkm false`. |

## Database / reference-path params

Most database paths (`--gtdbtk_db`, `--checkm2_db`, `--checkm1_db`,
`--singlem_metapackage`, `--sylph_db`, `--bakta_db`, `--dram_db`,
`--genomad_db`, `--checkv_db`, `--eggnog_db`, `--amrfinder_db`, `--host_ref`,
`--cleanifier_db`, `--cleanifier_nobjects`, `--genomespot_models`) are
documented with their tool and default in [databases.md](databases.md) — see
that table rather than duplicating it here.

Params not in that table:

| Param | Default | What it does |
|-------|---------|---------------|
| `--reference_genomes` | — | Directory of external reference genome FASTAs (metagenome modes). When set, references are dereplicated together with the HQ MAGs, reads are mapped to that combined set, and reference proteins are added to an expanded gene catalogue (+DRAM+RPKM). Requires `--skip_checkm false` and `--skip_dereplication false`. |
| `--reference_genome_extension` | `fasta,fa,fna,fasta.gz,fa.gz,fna.gz` | Comma-separated file extensions to glob under `--reference_genomes`. |
| `--reference_genomes_checkm2` | — | Optional CheckM2 `quality_report.tsv` covering the reference genomes; every reference must appear in it or the run fails. If unset, CheckM2 is run on the references directly. |
| `--download_db` | — | Comma-list of DB names to download (or `all`), written under `--db_outdir`. Only used with `--mode download_dbs`. |
| `--db_outdir` | `databases` | Where `--mode download_dbs` writes. |

## Tool-choice / algorithm params

| Param | Default | What it does |
|-------|---------|---------------|
| `--catalogue_identities` | `1.0,0.9` | Comma-separated CD-HIT identities for the gene catalogue; the first is the primary catalogue (used for CDS/membership/DRAM). |
| `--rpkm_min_read_length` | `140` | Minimum selected R1 read length retained for RPKM DIAMOND blastx. |
| `--sylph_profile_args` | — | Optional extra arguments for `sylph profile`. |
| `--marker_tree_builder` | `veryfasttree` | Tree builder for the metagenome marker-gene tree (`veryfasttree` or `iqtree`). |
| `--marker_tree_genome_source` | `representatives` | Which user genomes to place in the marker-gene tree: `representatives`, `hq_representatives` (HQ-first-then-dereplicated set, bypasses the completeness/contamination filter below), or `all_bins`. |
| `--marker_tree_min_completeness` | `90` | CheckM2 completeness threshold for genomes placed in the marker-gene tree. |
| `--marker_tree_max_contamination` | `5` | CheckM2 contamination threshold for genomes placed in the marker-gene tree. |
| `--marker_tree_use_closest` | `true` | Select closest reference genomes by GTDB-Tk tree topology. |
| `--marker_tree_use_related` | `true` | Select one reference per order sharing a bin's order but not its family. |
| `--marker_tree_closest_n` | `2` | Number of closest reference leaves per genome. |
| `--marker_tree_related_per_order` | `1` | Related references kept per order. |
| `--marker_tree_exclude_pattern` | — | Extra regex to exclude leaf names when selecting closest references. |
| `--marker_tree_reference_accessions` | — | Optional file of GB_/RS_ accessions to add to the marker-gene tree verbatim. |
| `--strain_genome_source` | `hq_representatives` | Which cross-sample genome set both strain tools map to: `hq_representatives` (HQ MAGs extracted from the full bin set then dereplicated), `hq_representatives_direct` (dereplicate first, then keep HQ representatives), `representatives` (any quality), or `hq_ref_representatives` (HQ MAGs dereplicated with external `--reference_genomes`). All `hq_*` options need a CheckM report. |
| `--strain_min_completeness` | `90` | CheckM completeness threshold for the strain reference set (MIMAG high-quality), applied on top of `--strain_genome_source`. Incompleteness biases calls toward "same strain", the conservative direction. Genomes absent from the CheckM report (external references) are kept. |
| `--strain_max_contamination` | `5` | CheckM contamination threshold for the strain reference set. Stricter than the pipeline HQ filter (`completeness - 3*contamination >= 50`) on purpose: contaminating contigs collect reads from unrelated populations and produce spurious SNPs, i.e. false "different strain" calls. |
| `--aviary_extra_binners` | — | Extra binner name(s) passed to Aviary recover's `--extra-binners` (e.g. `comebin` or `comebin maxbin2`). Metagenome modes only. |
| `--aviary_long_read_type` | `ont_hq` | Aviary `--long-read-type`. Only takes effect in `nanopore_metagenome` (no-op in `illumina_metagenome`, which always uses the paired-end read path). |
| `--long_read_type` | `ont_r10` | Autocycler read type. Used only by `nanopore_isolate`. |
| `--dorado_model` | `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` | Dorado basecalling model, used when basecalling POD5 input. |
| `--dorado_barcode_kit` | — | Optional Dorado barcode kit name for POD5 basecalling/demux. |
| `--dorado_device` | `auto` | Dorado `--device` value (`auto`, `cuda:0`, `cpu`, …). |
| `--comparison_manifest` | — | Optional manifest assigning sample/reference rows to comparison groups. See [usage.md](usage.md#isolate-comparison-manifest). Isolate modes only. |
| `--samples_include` | — | Optional sample id list used only for fallback samplesheet-group comparisons. Isolate modes only. |
| `--chewbbaca_training_file` | — | Optional chewBBACA/Prodigal training file. Isolate modes only. |
| `--chewbbaca_cgmlst_thresholds` | `0.90 0.95 0.99 1` | Space-separated thresholds passed to chewBBACA `ExtractCgMLST`. Isolate modes only. |

### Not currently implemented

These three are accepted (they validate and won't error) but have no effect
— the pipeline always uses the first listed tool regardless of the value set:

| Param | Enum | Currently always uses |
|-------|------|------------------------|
| `--catalogue_clusterer` | `cdhit` / `mmseqs` | CD-HIT (no MMseqs2 path exists in the gene-catalogue subworkflow) |
| `--core_alignment` | `parsnp` / `snippy` | ParSNP, gated by `--run_parsnp` (no Snippy path exists) |
| `--tree_builder` | `iqtree` / `raxml` | IQTree, gated by `--run_tree` (no RAxML path exists; not to be confused with `--marker_tree_builder`, which is wired up) |

## Resource params

| Param | Default | What it does |
|-------|---------|---------------|
| `--max_cpus` | `16` | CPU ceiling applied to every process via `process.resourceLimits`. |
| `--max_memory` | `128.GB` | Memory ceiling applied to every process. |
| `--max_time` | `72.h` | Time ceiling applied to every process. |

Finer-grained tuning (per-label CPU/memory in `conf/base.config`, or
tool-specific `ext.args` in `conf/modules.config` that aren't exposed as a
`--param`, e.g. dereplication ANI threshold, read-mapping identity
thresholds) isn't controlled by a flag — edit those config files directly.
See the "I want to change X" table in [configuration.md](configuration.md).

## Worked examples

```bash
# Metagenome run skipping SingleM/mobile-elements profiling, and disabling RPKM
nextflow run . -profile bunya --mode illumina_metagenome \
  --input samplesheet.csv --outdir results \
  --skip_singlem true --skip_mobile_elements true --skip_rpkm true

# Isolate run with only FastANI + chewBBACA in the comparative step
# (ParSNP and Panaroo, and the tree that depends on Panaroo, turned off)
nextflow run . -profile bunya --mode illumina_isolate \
  --input samplesheet.csv --outdir results \
  --run_parsnp false --run_panaroo false
```
