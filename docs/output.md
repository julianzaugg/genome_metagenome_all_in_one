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
pipeline_info/          # timeline / report / trace / dag
```

`high_quality_representatives/` holds bins passing completeness − 3×contamination ≥ 50
in **either** CheckM1 **or** CheckM2 (whichever ran). CoverM uses CheckM2 to pick
representatives, falling back to CheckM1 if CheckM2 is skipped.

`06_aviary/all_aviary_bins/` contains the sample-prefixed Aviary bins used by
CheckM, dereplication, read mapping, and GTDB-Tk. Filenames follow
`<sample>.<binner>.<bin_number>.fasta`, with `bin_contig_list.tsv` mapping each
renamed bin to its contigs.

The other modes get their own gapless `01..N` sequences when implemented (the
isolate-track numbers currently in `conf/modules.config` are placeholders for
the scaffolds).

## Provenance

- `09_dereplicated_bins/cluster_definition.tsv` — which bins collapsed into each
  representative.
- `12_gene_catalogue/gene_catalogue_membership.tsv` — which predicted gene
  (namespaced `<sample>___<gene>`) maps to each catalogue cluster.

Paths/numbers are set in `conf/modules.config` and easily changed.
