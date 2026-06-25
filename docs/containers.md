# Containers

This pipeline is **containers-only** (Apptainer). Two sources:

1. **Reused nf-core modules** (`modules/nf-core/*`) ship their own biocontainer
   directive — pulled from a registry and converted by Apptainer on first use.
   Override any of them with a `withName` block in `conf/containers.config`.

2. **Bespoke local modules** use local `.sif` images under `params.container_base`
   (set per profile: `/scratch/project/a_ace/containers` on Bunya,
   `/srv/db/containers` locally). Mapped in `conf/containers.config` as
   `"${params.container_base}/<name>.sif"`.

## Illumina-metagenome containers

### Auto-pulled from quay.io biocontainers — nothing to build
These are already wired in `conf/containers.config` to quay URIs; Apptainer pulls
and caches them on first use: **sylph, singlem, coverm, pyrodigal, cd-hit, dram,
checkm-genome (CheckM1), nonpareil, seqkit** (glue steps) and **python** (bin/
helper). Plus the nf-core modules (fastp, spades, gtdbtk, checkm2, geNomad, checkv).

### Must provide a local `.sif` at `params.container_base`
Either not on biocontainers, or a process needs two tools in one image:

| Image (`<name>.sif`)   | Why a local image | Used by |
|------------------------|-------------------|---------|
| `aviary_0.12.0`        | complex dependency stack | binning |
| `genomespot_1.0`       | not packaged on biocontainers | bin growth prediction (optional) |
| `minimap2_samtools`    | needs minimap2 **+** samtools together | host removal (optional) |
| `checkv_blast`         | needs blast+ **+** python3 together | virus/plasmid clustering |

For the two-tool images, either build a small combined image or use a mulled
biocontainer (e.g. nf-core's `minimap2`+`samtools` mulled image).

## Other-mode containers (scaffolds)
Filled as local `.sif` placeholders for now: `fastplong`, `cleanifier`, `myloasm`,
`autocycler`, `dorado`, `polypolish`, `dnaapler`, `shovill`, `chewbbaca`, `parsnp`,
`fastani`. Most (shovill, fastplong, polypolish, dnaapler, chewbbaca, parsnp,
fastani) are single biocontainers — swap them to quay URIs the same way when you
build those workflows. `dorado`, `myloasm`, `autocycler`, `cleanifier` are bespoke.

## Apptainer config

`apptainer.enabled = true` and `autoMounts = true` are set globally;
`apptainer.cacheDir` is set per profile. The `bunya_gpu` profile adds
`apptainer.runOptions = '--nv'` for GPU passthrough (dorado).
