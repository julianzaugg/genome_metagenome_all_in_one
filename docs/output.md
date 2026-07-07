# Output

Outputs are published under `--outdir` in numbered directories mirroring the
bash-pipeline convention. **Numbers run gaplessly in execution order, per mode** —
they are not kept consistent across modes (each `--mode` is a separate run with
its own output tree, so there's no reason to). Illumina metagenome layout:

```
00_read_stats/          # seqkit stats --tabular --all per stage + read_stat_report.tsv
01_fastp/               # QC'd reads + reports
02_sylph/               # sylph combined profile
03_singlem/             # multi-sample SingleM profile + OTU table
04_host_removed/        # host-filtered reads
05_metaspades/          # assemblies (scaffolds)
06_aviary/              # Aviary recovery; all_aviary_bins/ is the renamed canonical bin set
07_checkm2/             # CheckM2 on all bins (drives dereplication + HQ selection)
08_dereplicated_bins/   # CoverM cluster: representatives/ + high_quality_representatives/ + cluster_definition.tsv
08_dereplicated_hq_bins/ # CoverM cluster on HQ-first bins: HQ MAGs extracted from the FULL bin set, THEN dereplicated
08_dereplicated_hq_ref_bins/ # HQ MAGs dereplicated TOGETHER with external reference genomes (if --reference_genomes)
09_coverm_bins/         # per-sample abundance + bam vs all dereplicated representatives
09_coverm_hq_bins/      # per-sample abundance + bam vs high_quality_representatives/ only (no competing siblings)
09_coverm_hq_derep_bins/ # per-sample abundance + bam vs the HQ-first-then-dereplicated set (08_dereplicated_hq_bins)
09_coverm_hq_ref_bins/  # per-sample abundance + bam vs the HQ-MAGs+references set (08_dereplicated_hq_ref_bins; if --reference_genomes)
10_coverm_scaffolds/    # per-sample coverage/counts + bam vs assembled scaffolds
11_pyrodigal/           # predicted proteins/genes per assembly
12_gene_catalogue/      # cd-hit catalogue(s) + nucleotide CDS + membership (provenance)
12_gene_catalogue_expanded/ # as 12_gene_catalogue but scaffold + reference-genome proteins (if --reference_genomes)
13_dram/                # DRAM functional annotation of the catalogue (annotated in parallel chunks, merged)
13_dram_expanded/       # DRAM functional annotation of the expanded catalogue (if --reference_genomes)
14_dram_bins/           # per-bin DRAM annotation + combined cross-MAG distillate in all_bins/ (if --run_dram_bins)
15_gtdbtk/              # GTDB-Tk classification of representatives
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
25_reference_genomes/   # normalised reference FASTAs, their CheckM2 report, predicted proteins (if --reference_genomes)
pipeline_info/          # timeline / report / trace / dag
```

`24_marker_tree/` (opt-in via `--run_marker_tree`) holds, per domain
(`bac120`/`ar53`): `<domain>.marker_msa.fasta` (placed genomes + selected
references, from the GTDB-Tk marker alignment), `<domain>.treefile` (VeryFastTree
by default, or IQ-TREE via `--marker_tree_builder iqtree`), and the
`<domain>.reference_genomes.tsv` / `<domain>.closest_references.tsv` selection
tables. Placed genomes are the `--marker_tree_genome_source` set filtered by
CheckM2 `--marker_tree_min_completeness` / `--marker_tree_max_contamination`.
References come from closest-by-topology (`--marker_tree_use_closest`) and/or
same-order-different-family (`--marker_tree_use_related`) selection, plus any
`--marker_tree_reference_accessions`; disable both selection modes for a
genomes-only tree.

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
  (`08_dereplicated_hq_ref_bins/`).

  `09_coverm_hq_derep_bins/` is a third HQ mapping against a **differently
  constructed** set (`08_dereplicated_hq_bins/`). The `_HQ_MAGs*` sets above
  dereplicate the full bin set first and then keep the representatives that pass
  the HQ filter, so a cluster whose chosen representative isn't HQ contributes no
  HQ MAG. `08_dereplicated_hq_bins/` reverses the order — HQ MAGs are extracted
  from the FULL pre-dereplication bin set first, then those are dereplicated — so
  every HQ cluster is represented by an HQ genome. Which CheckM report(s) decide
  "HQ" is controlled by `--hq_quality_source` (`both` (default) | `checkm1` |
  `checkm2`; a bin is HQ if it passes in any selected report).
- **Isolate**: same QC-stage columns, then CoverM mapping stats against the
  sample's own assembly — `Covered_fraction, Mean_coverage, Read_count,
  Read_count_percent` (with `_SR` variants when a hybrid nanopore isolate also
  has short reads mapped).

`GBbp` is the raw total bases (forward + reverse) in gigabasepairs.

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
