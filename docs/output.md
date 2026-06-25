# Output

Outputs are published under `--outdir` in numbered directories mirroring the
bash-pipeline convention. **Numbers run gaplessly in execution order, per mode** —
they are not kept consistent across modes (each `--mode` is a separate run with
its own output tree, so there's no reason to). Illumina metagenome layout:

```
01_fastp/               # QC'd reads + reports
02_sylph/               # sylph combined profile
03_singlem/             # per-sample singlem profiles + OTU tables
04_host_removed/        # host-filtered reads
05_metaspades/          # assemblies (scaffolds)
06_aviary/              # Aviary recovery (bins)
07_checkm2/             # bin completeness/contamination (drives dereplication)
08_dereplicated_bins/   # CoverM cluster: representatives/ + cluster_definition.tsv
09_coverm_bins/         # per-sample abundance vs representatives
10_pyrodigal/           # predicted proteins/genes per assembly
11_gene_catalogue/      # cd-hit catalogue + membership table (provenance)
12_dram/                # DRAM functional annotation of the catalogue
13_gtdbtk/              # GTDB-Tk classification of representatives
14_checkm1/             # CheckM1 (if --run_checkm1)
15_nonpareil/           # coverage redundancy (if --run_nonpareil)
16_genomespot/          # growth predictions (if --run_genomespot)
17_genomad/             # virus/plasmid prediction + pooled sequences
18_checkv/              # CheckV quality of pooled viruses
19_checkv_clustering/   # ANI clusters (virus + plasmid)
pipeline_info/          # timeline / report / trace / dag
```

The other modes get their own gapless `01..N` sequences when implemented (the
isolate-track numbers currently in `conf/modules.config` are placeholders for
the scaffolds).

## Provenance

- `09_dereplicated_bins/cluster_definition.tsv` — which bins collapsed into each
  representative.
- `12_gene_catalogue/gene_catalogue_membership.tsv` — which predicted gene
  (namespaced `<sample>___<gene>`) maps to each catalogue cluster.

Paths/numbers are set in `conf/modules.config` and easily changed.
