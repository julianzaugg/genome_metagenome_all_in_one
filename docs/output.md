# Output

Outputs are published under `--outdir` in numbered directories mirroring the
bash-pipeline convention. **Numbers run gaplessly in execution order, per mode** —
they are not kept consistent across modes (each `--mode` is a separate run with
its own output tree, so there's no reason to). Illumina metagenome layout:

```
00_read_stats/          # seqkit stats --tabular --all for raw and QC-step reads
01_fastp/               # QC'd reads + reports
02_sylph/               # sylph combined profile
03_singlem/             # multi-sample SingleM profile + OTU table
04_host_removed/        # host-filtered reads
05_metaspades/          # assemblies (scaffolds)
06_aviary/              # Aviary recovery; all_aviary_bins/ is the renamed canonical bin set
07_checkm2/             # CheckM2 on all bins (drives dereplication + HQ selection)
08_dereplicated_bins/   # CoverM cluster: representatives/ + high_quality_representatives/ + cluster_definition.tsv
09_coverm_bins/         # per-sample abundance vs representatives
10_coverm_scaffolds/    # per-sample coverage/counts vs assembled scaffolds
11_pyrodigal/           # predicted proteins/genes per assembly
12_gene_catalogue/      # cd-hit catalogue(s) + nucleotide CDS + membership (provenance)
13_dram/                # DRAM functional annotation of the catalogue
14_gtdbtk/              # GTDB-Tk classification of representatives
15_checkm1/             # CheckM1 on all bins (if --run_checkm1) — also feeds HQ selection
16_nonpareil/           # coverage redundancy (if --run_nonpareil)
17_genomespot/          # growth predictions (if --run_genomespot)
18_barrnap/             # rRNA / 16S per representative (if --run_barrnap)
19_genomad/             # virus/plasmid prediction + pooled seqs/proteins/genes/summary
20_checkv/              # CheckV quality of pooled viruses
21_checkv_clustering/   # ANI clusters (virus + plasmid)
22_rpkm/                # SingleM-normalized RPKM for the gene catalogue
pipeline_info/          # timeline / report / trace / dag
```

`high_quality_representatives/` holds bins passing completeness − 3×contamination ≥ 50
in **either** CheckM1 **or** CheckM2 (whichever ran). CoverM uses CheckM2 to pick
representatives, falling back to CheckM1 if CheckM2 is skipped.

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
09_coverm_bins/         # minimap2-ont, 90% identity
10_coverm_scaffolds_nanopore/
```

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
14_panaroo/
15_parsnp/
16_gubbins/
17_fastani/
18_chewbacca/
19_tree/
```

Nanopore isolate mapping emits separate CoverM outputs for long reads
(`10_coverm_scaffolds_nanopore/`) and hybrid Illumina reads
(`10_coverm_scaffolds_illumina/`) when both are present.

## Provenance

- `09_dereplicated_bins/cluster_definition.tsv` — which bins collapsed into each
  representative.
- `12_gene_catalogue/gene_catalogue_membership.tsv` — which predicted gene
  (namespaced `<sample>___<gene>`) maps to each catalogue cluster.
- `22_rpkm/gene_catalogue_rpkm_per_gene_normalised.tsv` — gene-catalogue RPKM
  normalized by each sample's mean SingleM marker-gene RPKM.
- `22_rpkm/gene_catalogue_mapped_reads_per_gene.tsv` — selected R1 DIAMOND read
  counts per catalogue gene and sample.
- `22_rpkm/singlem_sample_rpkm.tsv` and `22_rpkm/singlem_rpkm_means.tsv` —
  marker-level and per-sample SingleM normalization values.
- `14_comparison_groups/<group>/entries.tsv` — isolate comparison membership,
  including references and the chosen Parsnp reference.
- `18_chewbacca/<group>_chewbbaca/genome_hash_map.tsv` — mapping from original
  genome ids to chewBBACA-safe FASTA names.

Paths/numbers are set in `conf/modules.config` and easily changed.
