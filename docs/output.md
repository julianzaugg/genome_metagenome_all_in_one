# Output

Outputs are published under `--outdir` in numbered directories mirroring the
bash-pipeline convention. **Numbers run gaplessly in execution order, per mode** —
they are not kept consistent across modes (each `--mode` is a separate run with
its own output tree, so there's no reason to). Illumina metagenome layout:

```
00_read_stats/          # seqkit stats --tabular --all per stage + read_stat_report.tsv + mapping_assessment.tsv (long-format)
01_fastp/               # QC'd reads + reports
02_sylph/               # sylph combined profile
03_singlem/             # multi-sample SingleM profile + OTU table
04_host_removed/        # host-filtered reads
05_metaspades/          # assemblies (scaffolds)
06_aviary/              # Aviary recovery; all_aviary_bins/ is the renamed canonical bin set
07_checkm2/             # CheckM2 on all bins (drives dereplication + HQ selection)
08_dereplicated_bins/   # CoverM cluster: representatives/ + high_quality_representatives/ + cluster_definition.tsv
08_dereplicated_hq_bins/ # CoverM cluster on HQ-first bins: HQ MAGs extracted from the FULL bin set, THEN dereplicated
08_dereplicated_hq_ref_bins/ # HQ MAGs (extracted from the FULL bin set, independently of 08_dereplicated_hq_bins) dereplicated TOGETHER with external reference genomes (if --reference_genomes)
08_within_sample_dereplicated_bins/<id>/    # per-sample(or -group) dereplicated representatives (if --within_sample_dereplication sample|group)
08_within_sample_dereplicated_hq_bins/<id>/ # per-sample(or -group) HQ-first-then-dereplicated MAGs
09_coverm_bins/         # per-sample abundance + bam + mapping_assessment.tsv[_per_genome] vs all dereplicated representatives
09_coverm_hq_bins/      # per-sample abundance + bam + mapping_assessment.tsv[_per_genome] vs high_quality_representatives/ only (no competing siblings)
09_coverm_hq_derep_bins/ # per-sample abundance + bam + mapping_assessment.tsv[_per_genome] vs the HQ-first-then-dereplicated set (08_dereplicated_hq_bins)
09_coverm_hq_ref_bins/  # per-sample abundance + bam + mapping_assessment.tsv[_per_genome] vs the HQ-MAGs+references set (08_dereplicated_hq_ref_bins; if --reference_genomes)
09_coverm_within_sample_derep_bins/ # each sample's reads + mapping_assessment.tsv[_per_genome] vs its own within-sample(or -group) dereplicated bins
09_coverm_within_sample_hq_bins/     # each sample's reads + mapping_assessment.tsv[_per_genome] vs its own within-sample(or -group) HQ MAGs
10_coverm_scaffolds/    # per-sample coverage/counts + bam + mapping_assessment.tsv[_per_genome] vs assembled scaffolds
11_pyrodigal/           # predicted proteins/genes per assembly
12_gene_catalogue/      # cd-hit catalogue(s) + nucleotide CDS + membership (provenance)
12_gene_catalogue_expanded/ # as 12_gene_catalogue but scaffold + reference-genome proteins (if --reference_genomes)
13_dram/                # DRAM functional annotation of the catalogue (annotated in parallel chunks, merged)
13_dram_expanded/       # DRAM functional annotation of the expanded catalogue (if --reference_genomes)
14_dram_bins/           # per-bin DRAM annotation + combined cross-MAG distillate in all_bins/ (if --run_dram_bins)
15_gtdbtk/              # GTDB-Tk classification of all bins (+ references if --reference_genomes_taxonomy / --marker_tree_include_references)
16_checkm1/             # CheckM1 on all bins (if --run_checkm1) — also feeds HQ selection
17_nonpareil/           # coverage redundancy (if --run_nonpareil)
18_genomespot/          # growth predictions (if --run_genomespot)
19_barrnap/             # rRNA / 16S per representative (if --run_barrnap)
20_genomad/             # virus/plasmid prediction + pooled seqs/proteins/genes/summary
21_checkv/              # CheckV quality of pooled viruses
22_checkv_clustering/   # ANI clusters (virus + plasmid)
23_rpkm/                # SingleM-normalized RPKM for the gene catalogue
23_rpkm_expanded/       # RPKM for the expanded catalogue; reuses 23_rpkm's SingleM marker blast (if --reference_genomes)
24_marker_tree/         # MAG + GTDB-reference marker-gene tree (if --run_marker_tree)
25_reference_genomes/   # normalised reference FASTAs, their CheckM2 report, predicted proteins, USERREF_-prefixed copies for GTDB-Tk (if --reference_genomes)
26_strain_reference/    # genomes used for strain comparison + audit table, combined FASTA, .stb (if --run_instrain / --run_tracs)
27_instrain/            # inStrain profiles, compare output, and strain-sharing summary tables (if --run_instrain)
28_tracs/               # TRACS reference db, pairwise SNP/transmission distances, strain clusters (if --run_tracs)
pipeline_info/          # timeline / report / trace / dag
```

`24_marker_tree/` (opt-in via `--run_marker_tree`) holds, per domain
(`bac120`/`ar53`): `<domain>.marker_msa.fasta` (placed genomes + selected
references, from the GTDB-Tk marker alignment), `<domain>.treefile` (VeryFastTree
by default, or IQ-TREE via `--marker_tree_builder iqtree`), and the
`<domain>.reference_genomes.tsv` / `<domain>.closest_references.tsv` selection
tables. Placed genomes are the `--marker_tree_genome_source` set
(`representatives` — the default, `hq_representatives` — the HQ-first-then-
dereplicated set (HQ bins, `completeness − 3×contamination ≥ 50`, extracted
from the full bin set and clustered among themselves), or `all_bins`) filtered by CheckM2
`--marker_tree_min_completeness` / `--marker_tree_max_contamination` (the numeric
filter is skipped for `hq_representatives`, which is already high-quality). GTDB
*database* references come from closest-by-topology (`--marker_tree_use_closest`)
and/or same-order-different-family (`--marker_tree_use_related`) selection, plus any
`--marker_tree_reference_accessions`; disable both selection modes for a
genomes-only tree. With `--marker_tree_include_references`, the external
`--reference_genomes` are classified by GTDB-Tk (renamed `USERREF_*` internally to
avoid GCA_/GCF_ accession collisions) and placed in the tree alongside the MAGs,
regardless of their CheckM2 quality; their leaf labels and the published GTDB-Tk
summary show the original names. `--reference_genomes_taxonomy` classifies the
references without placing them in a tree.

`high_quality_representatives/` holds bins passing completeness − 3×contamination ≥ 50
in **either** CheckM1 **or** CheckM2 (whichever ran). CoverM uses CheckM2 to pick
representatives, falling back to CheckM1 if CheckM2 is skipped.

`00_read_stats/read_stat_report.tsv` tracks per-sample read counts through every
step that ran. Counts are summed across mate files (paired reads = forward +
reverse), and **every `*_percent` column is relative to the raw input**. Column
shape depends on the mode:

- **Metagenome**: `Sample_ID, GBbp, Raw_count`, then a `<Tool>_count`/`<Tool>_percent`
  pair per QC stage (Fastp / Porechop / Fastplong / Cleanifier — whichever ran),
  then `Reads_mapped_Scaffolds`, `Reads_mapped_Dereplicated_Bins`,
  `Reads_mapped_HQ_MAGs`, `Reads_mapped_HQ_MAGs_direct`, and
  `Reads_mapped_HQ_Derep_MAGs` count/percent (from the CoverM `Count` method).
  `Dereplicated_Bins` is the full set of cluster representatives of any quality;
  `MAGs` is reserved for the high-quality subset. The HQ columns can differ:
  `_HQ_MAGs` is the subset of the full `09_coverm_bins/` mapping whose genome
  matches `high_quality_representatives/`, while `_HQ_MAGs_direct` comes from a
  separate `09_coverm_hq_bins/` mapping against the HQ MAGs only. When
  dereplication leaves redundant near-identical lower-quality bins alongside an
  HQ rep, competitive mapping in the full-set run splits reads across those
  siblings and the subset count undercounts — the direct mapping doesn't have
  that problem. `_HQ_Derep_MAGs` is the `09_coverm_hq_derep_bins/` mapping
  against the HQ-first-then-dereplicated set (`08_dereplicated_hq_bins/`; see
  below). With `--reference_genomes`, a `Reads_mapped_HQ_Ref_MAGs` count/percent
  pair is also added — the `09_coverm_hq_ref_bins/` mapping against the HQ MAGs
  dereplicated together with the external reference genomes
  (`08_dereplicated_hq_ref_bins/`). This HQ extraction is its own independent
  pass over the full pre-dereplication bin set — it does not reuse the HQ MAGs
  from `08_dereplicated_hq_bins/` — so it re-applies the same HQ filter before
  adding the references and clustering.

  `09_coverm_hq_derep_bins/` is a third HQ mapping against a **differently
  constructed** set (`08_dereplicated_hq_bins/`). The `_HQ_MAGs*` sets above
  dereplicate the full bin set first and then keep the representatives that pass
  the HQ filter, so a cluster whose chosen representative isn't HQ contributes no
  HQ MAG. `08_dereplicated_hq_bins/` reverses the order — HQ MAGs are extracted
  from the FULL pre-dereplication bin set first, then those are dereplicated — so
  every HQ cluster is represented by an HQ genome. Which CheckM report(s) decide
  "HQ" is controlled by `--hq_quality_source` (`both` (default) | `checkm1` |
  `checkm2`; a bin is HQ if it passes in any selected report).

  With `--within_sample_dereplication sample` (or `group`), two more pairs appear:
  `Reads_mapped_PerSample_Derep_MAGs` and `Reads_mapped_PerSample_HQ_MAGs`. These
  map each sample's reads to bins dereplicated **within that sample (or group)
  only** (`09_coverm_within_sample_derep_bins/` and `09_coverm_within_sample_hq_bins/`),
  rather than across all samples. This scope is independent of `--skip_dereplication`
  (across-all): run either alone, or both together.
- **Isolate**: same QC-stage columns, then CoverM mapping stats against the
  sample's own assembly — `Covered_fraction, Mean_coverage, Read_count,
  Read_count_percent` (with `_SR` variants when a hybrid nanopore isolate also
  has short reads mapped).

`GBbp` is the raw total bases (forward + reverse) in gigabasepairs.

### Bases-mapped metrics (unless `--skip_mapping_assessment`)

A read count is a fine proxy for "how much data is represented" only when every
read is roughly the same length. For nanopore reads, which span three orders of
magnitude in length, it is not: the *number* of reads that mapped says little
about the *fraction of sequenced bases* a genome set accounts for. Alongside
every `Reads_mapped_<Label>_count/percent` pair above, the report also carries
three base-level metrics, computed by `MAPPING_ASSESS`
(`bin/assess_mapping.py`) from the same CoverM-cached BAM, re-filtered with
`coverm filter` at the same identity/aligned-percent thresholds so these
numbers describe the same alignment set as the `Count` column:

| Column | Metric | Definition |
| --- | --- | --- |
| `Bases_mapped_A_<Label>_count/percent` | **A** | Full SEQ length of every primary mapped read (a read counted once at its full length — "how much read data mapped somewhere"). The number to lead with for long reads. |
| `Bases_mapped_B_<Label>_count/percent` | **B** | Aligned bases (CIGAR `M`/`I`/`=`/`X`) of **primary alignments only** — a stricter, lower-bound figure that excludes soft-clipped/unaligned tails. |
| `Bases_mapped_C_<Label>_count/percent` | **C** | Aligned bases of primary **and supplementary** alignments — the honest total for a long read that splits across a contig/genome boundary or is chimeric. A read whose primary lands on one genome but has a supplementary segment aligning to another contributes bases to both genomes' per-genome tables (see below), even though `Count` credits the whole read to the primary's genome only. |

All three divide by the **raw sequenced bases** (the same total `GBbp` is
derived from) — a different denominator from the `_percent` columns next to
`Count`, which divide by raw *reads*. `Bases_mapped_*_HQ_MAGs` (the subset
column, alongside `Reads_mapped_HQ_MAGs`) is summed from the `Dereplicated_Bins`
per-genome table over genomes in `high_quality_representatives/`, mirroring how
`Reads_mapped_HQ_MAGs` is a row-subset sum of the `09_coverm_bins/` table.

Two more files land alongside each mapping's existing abundance table:
`{sample}.{Label}.mapping_assessment.tsv` (the one-row source of the columns
above) and `{sample}.{Label}.mapping_assessment_per_genome.tsv` (the same three
metrics broken down per genome/contig — the source of the `HQ_MAGs` subset
sum, and useful on its own for seeing which genomes actually received the
bases from a split long read). `00_read_stats/mapping_assessment.tsv`
concatenates every sample's `.mapping_assessment.tsv` into one long-format
table (`sample, set, ...`) across all reference sets that ran.

**Caveat:** `--min-read-aligned-percent 0.75` on every ONT read-mapping step
(`conf/modules.config`) is evaluated per **primary** alignment record
(aligned bases / full read length). A chimeric or split long read whose
primary segment alone covers less than 75% of the read is discarded entirely
by CoverM — before `Count`, and before metric C ever sees it — which is
exactly the kind of read metric C exists to give credit to. This threshold is
left unchanged for result comparability; lowering or removing it (and
rerunning) would let both `Count` and metric C credit these reads. See the
comment beside the `minimap2-ont` `ext.args` entries in `conf/modules.config`.

`06_aviary/all_aviary_bins/` contains the sample-prefixed Aviary bins used by
CheckM, dereplication, read mapping, and GTDB-Tk. Filenames follow
`<sample>.<binner>.<bin_number>.fasta`, with `bin_contig_list.tsv` mapping each
renamed bin to its contigs.

Nanopore metagenome reuses the metagenome output families, with Nanopore-specific
read QC and assembly directories:

```
01_dorado_basecall/      # only when POD5 input is basecalled
02_porechop/
03_fastplong/
07_myloasm/
09_coverm_bins/         # minimap2-ont, 90% identity, + bam
09_coverm_hq_bins/      # minimap2-ont, 90% identity, HQ representatives only, + bam
09_coverm_hq_derep_bins/ # minimap2-ont, HQ-first-then-dereplicated set, + bam
09_coverm_hq_ref_bins/  # minimap2-ont, HQ-MAGs+references set, + bam (if --reference_genomes)
09_coverm_within_sample_derep_bins/ # minimap2-ont, within-sample(or -group) dereplicated set, + bam
09_coverm_within_sample_hq_bins/     # minimap2-ont, within-sample(or -group) HQ MAGs, + bam
10_coverm_scaffolds_nanopore/ # + bam
```

The `--reference_genomes` feature applies to both metagenome tracks: references are
dereplicated with the HQ MAGs (`08_dereplicated_hq_ref_bins/`), reads are mapped to
that set (`09_coverm_hq_ref_bins/`), and reference proteins are added to an expanded
gene catalogue (`12_gene_catalogue_expanded/` + `13_dram_expanded/`). RPKM for the
expanded catalogue (`23_rpkm_expanded/`) is Illumina-only, and reuses the SingleM
marker blast from `23_rpkm/` (only the gene-catalogue blast is recomputed). References
are scored with CheckM2 (a report is generated, or supply one with
`--reference_genomes_checkm2`, in which case every reference must appear in it).

Isolate workflows publish assembly, QC, annotation, mobile-element, mapping, and
comparative outputs under their tool names. Key directories:

```
05_autocycler/           # Nanopore isolate assembly
06_dorado_polish/        # Nanopore isolate polishing, if enabled
07_shovill/              # Illumina isolate assembly
07_polypolish/           # hybrid Nanopore isolate polishing, if short reads exist
08_dnaapler/             # Nanopore isolate chromosome orientation
10_bakta/
11_mlst/
12_amrfinder/
13_isescan/
14_comparison_groups/    # materialized sample/reference groups
15_panaroo/
16_parsnp/
17_gubbins/
18_fastani/
19_chewbacca/
20_tree/
```

Nanopore isolate mapping emits separate CoverM outputs for long reads
(`10_coverm_scaffolds_nanopore/`) and hybrid Illumina reads
(`10_coverm_scaffolds_illumina/`) when both are present.


## Strain comparison (`26_strain_reference/`, `27_instrain/`, `28_tracs/`)

Opt-in via `--run_instrain` (illumina_metagenome only) and/or `--run_tracs` (both
metagenome modes). These answer a question dereplication cannot: 95% ANI clustering
establishes that two samples share a **species**; these tools establish whether they share
a **strain**.

Every sample's reads are mapped to **one shared cross-sample reference set**, and the
resulting per-sample allele profiles are compared pairwise. A shared reference is the whole
point — per-sample references would give each sample its own coordinate system and nothing
would be comparable. This is why both tools require `--skip_dereplication false`.
`--within_sample_dereplication` is a separate, additive path (it clusters each sample's bins
on their own) and is neither a substitute nor an obstacle: it can be on or off.

`26_strain_reference/` holds the shared reference. `strain_reference_genomes.tsv` is the
audit table — one row per candidate genome with its completeness, contamination, which
CheckM report scored it, and whether it was kept. The set starts from
`--strain_genome_source` (default `hq_representatives`, the HQ-first-then-dereplicated set)
and is then filtered to `--strain_min_completeness` / `--strain_max_contamination`
(default 90/5, MIMAG high-quality). That is deliberately stricter than the pipeline's HQ
filter (`completeness − 3×contamination ≥ 50`), which admits e.g. a 95%-complete /
15%-contaminated bin — fine for abundance, too loose here (see the caveats below). Genomes
absent from the CheckM report — external `--reference_genomes` — are kept. If this table
shows almost everything dropped, that is a real signal about MAG quality, not a bug; the two
threshold params are the dial. `strain_reference.fasta` / `.stb` are the combined reference
for inStrain, with every contig header prefixed by its bin name (bins from separate
per-sample assemblies can otherwise share `NODE_..._length_..._cov_...` names).

`27_instrain/` (inStrain, short reads only). `profiles/<sample>.IS/` are the per-sample
microdiversity profiles; `strain_compare.IS/output/` holds the raw pairwise tables. The key
statistic is **popANI**, which counts a position as a difference only when the two samples
share *no* allele there — so a site where one sample is fixed and the other is polymorphic
for the same base is not a difference. That is what makes it a strain-sharing statistic
rather than a consensus-similarity one. It is only meaningful over a decent
`percent_genome_compared`, hence the conventional call: popANI ≥ 0.99999 over ≥ 50% of the
genome (set by `-ani`/`-cov` in `conf/modules.config`).

`summary/` reshapes that into the tables to actually read:
- `strain_sharing_counts.tsv` — per sample pair, how many genomes were comparable and how
  many shared a strain. Start here.
- `strain_sharing_summary.tsv` — one row per (genome, sample pair) with popANI, conANI,
  percent_genome_compared, cluster ids and the same/different call.
- `strain_sharing_matrix/<genome>.popani.tsv` — sample × sample popANI, for plotting.

The same/different call is taken from inStrain's own `strain_clusters` output rather than
re-derived, so tuning `-ani`/`-cov` changes the call in exactly one place.

`28_tracs/` (TRACS, short **or** long reads). `transmission_distances.csv` gives pairwise
SNP distances per reference genome (with a `filtered SNP distance` column and a
`sites considered` count), and `strain_clusters.csv` groups samples into transmission
clusters by single linkage. TRACS uses an empirical Bayes model over variable coverage
rather than fixed thresholds, and its distances are explicit **lower bounds**. Its reference
database (`strain_db.zip`) is built from the pipeline's own MAGs via `tracs build-db`, which
embeds both the sourmash index and the genomes themselves — so no GTDB download and no
network access at run time.

**inStrain is Illumina-only by design.** Its SNV model assumes short, low-error reads, and
there is no validated long-read parameterisation — swapping in a long-read aligner makes it
run without making the statistics sound. `--run_instrain` errors in `nanopore_metagenome`;
use `--run_tracs`, which supports `map-ont` natively.

### Reading these results honestly

- **Depth is the binding constraint, not the reference.** A genome needs adequate coverage
  in *both* samples of a pair to be compared at all. A missing row means "not enough data",
  which is easy to misread as "different strains".
- **Only species you binned are visible.** A species present in two samples but assembled in
  neither never appears. This is the cost of using your own MAGs instead of a full GTDB
  reference set — and it applies identically to both tools, so their calls stay comparable.
- **The two failure modes push in opposite directions.** Incompleteness shrinks the
  comparable fraction of a genome, so you see proportionally fewer differences and drift
  toward "same strain" — conservative, and visible in `percent_genome_compared` /
  `sites considered`. Contamination is the dangerous one: foreign contigs collect reads from
  unrelated populations, and when those populations differ between samples they generate
  spurious SNPs and drift toward "different strain". That asymmetry is why the reference set
  is gated on contamination directly, and why `TRACS_DISTANCE` keeps `--filter` on.
- **Mild asymmetric reference bias is inherent.** If a MAG was assembled from sample A, A's
  reads map to their own assembly near-perfectly while B's map to a foreign reference.
  inStrain's `--database_mode` exists to acknowledge exactly this. It is standard for any
  map-to-dereplicated-MAGs design, not a defect of this setup.

## Provenance

- `08_dereplicated_bins/cluster_definition.tsv` — which bins collapsed into each
  representative.
- `12_gene_catalogue/gene_catalogue_membership.tsv` — which predicted gene
  (namespaced `<sample>___<gene>`) maps to each catalogue cluster.
- `23_rpkm/gene_catalogue_rpkm_per_gene_normalised.tsv` — gene-catalogue RPKM
  normalized by each sample's mean SingleM marker-gene RPKM.
- `23_rpkm/gene_catalogue_mapped_reads_per_gene.tsv` — selected R1 DIAMOND read
  counts per catalogue gene and sample.
- `23_rpkm/singlem_sample_rpkm.tsv` and `23_rpkm/singlem_rpkm_means.tsv` —
  marker-level and per-sample SingleM normalization values.
- `14_comparison_groups/<group>/entries.tsv` — isolate comparison membership,
  including references and the chosen Parsnp reference.
- `19_chewbacca/<group>_chewbbaca/genome_hash_map.tsv` — mapping from original
  genome ids to chewBBACA-safe FASTA names.

Paths/numbers are set in `conf/modules.config` and easily changed.
