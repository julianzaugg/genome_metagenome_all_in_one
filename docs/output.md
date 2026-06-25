# Output

Outputs are published under `--outdir` in numbered directories mirroring the
bash-pipeline convention. Illumina metagenome layout:

```
03_fastp/               # QC'd reads + reports
04_sylph/               # sylph combined profile
05_singlem/             # per-sample singlem profiles + OTU tables
06_host_removed/        # host-filtered reads
07_metaspades/          # assemblies (scaffolds)
08_aviary/              # Aviary recovery (bins)
09_dereplicated_bins/   # CoverM cluster: representatives/ + cluster_definition.tsv
10_coverm_bins/         # per-sample abundance vs representatives
11_pyrodigal/           # predicted proteins/genes per assembly
12_gene_catalogue/      # cd-hit catalogue + membership table (provenance)
13_dram/                # DRAM functional annotation of the catalogue
14_gtdbtk/              # GTDB-Tk classification of representatives
15_checkm2/  15_checkm1/ # bin completeness/contamination
16_nonpareil/           # coverage redundancy (if --run_nonpareil)
17_genomespot/          # growth predictions (if --run_genomespot)
18_genomad/             # virus/plasmid prediction + pooled sequences
19_checkv/              # CheckV quality of pooled viruses
20_checkv_clustering/   # ANI clusters (virus + plasmid)
pipeline_info/          # timeline / report / trace / dag
```

Isolate modes reuse overlapping numbers with isolate-specific steps
(`10_bakta`, `11_mlst`, `12_amrfinder`, `13_isescan`, `14_panaroo`,
`15_parsnp`, `16_gubbins`, `17_fastani`, `18_chewbacca`, `19_tree`).

## Provenance

- `09_dereplicated_bins/cluster_definition.tsv` — which bins collapsed into each
  representative.
- `12_gene_catalogue/gene_catalogue_membership.tsv` — which predicted gene
  (namespaced `<sample>___<gene>`) maps to each catalogue cluster.

Paths/numbers are set in `conf/modules.config` and easily changed.
