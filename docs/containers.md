# Containers

This pipeline is **containers-only** (Apptainer). Two sources:

1. **Reused nf-core modules** (`modules/nf-core/*`) ship their own biocontainer
   directive — pulled from a registry and converted by Apptainer on first use.
   Override any of them with a `withName` block in `conf/containers.config`.

2. **Bespoke local modules** use local `.sif` images under `params.container_base`.
   The `local` and `bunya` profiles default this to `${projectDir}/containers`.
   Mapped in `conf/containers.config` as `"${params.container_base}/<name>.sif"`.

### Auto-pulled from quay.io biocontainers — nothing to build
Nearly every tool is wired to a quay.io biocontainer in `conf/containers.config`
(exact build tags resolved from the quay API); Apptainer pulls and caches each on
first use. This covers — across all modes — sylph, singlem, fastplong, cleanifier,
shovill, myloasm, autocycler, polypolish, dnaapler, **aviary**, coverm, pyrodigal,
cd-hit, mmseqs2, dram, checkm-genome (CheckM1), nonpareil, blast (checkv
clustering), bakta, chewbbaca, parsnp, fastani, seqkit (glue), python (helpers),
and the DB downloaders. Plus the nf-core modules (fastp, spades, gtdbtk, checkm2,
geNomad, checkv).

### Must provide a local `.sif` at `params.container_base`
Only tools not packaged on biocontainers need local images:

| Image (`<name>.sif`)   | Why a local image | Used by |
|------------------------|-------------------|---------|
| `dorado_1.4.0`         | ONT-proprietary, not on biocontainers | Nanopore basecall/polish |
| `genomespot_1.0`       | not packaged on biocontainers | bin growth prediction (optional) |

For the **Illumina-metagenome path you need no local builds** — host removal uses
the `cleanifier` biocontainer. Supply either a prebuilt `--cleanifier_db` index
or a FASTA with `--host_ref` so the pipeline can build the index. For Nanopore
you'll need `dorado`.

> These local-`.sif` entries use `${params.container_base}`, which is only set by a
> profile (e.g. `-profile local` → `<projectDir>/containers`). If you run without
> that profile, `container_base` is null and the path resolves to `null/...` — so
> always select the profile, or set `--container_base` explicitly, when a step
> needs one.

Caveats:
- **Aviary** biocontainer carries the CLI but builds its own tool conda envs at
  runtime — on an offline node you may prefer a self-contained image; override the
  `AVIARY_.*` line if so.
- **CHECKV_CLUSTER** uses the `blast` biocontainer plus the vendored stdlib
  `anicalc.py`/`aniclust.py`; if that image lacks `python3`, point it at a small
  blast+python image instead.

## Apptainer config

`apptainer.enabled = true` and `autoMounts = true` are set globally. The `local`
and `bunya` profiles set `apptainer.cacheDir` to `${projectDir}/.apptainer_cache`,
so auto-pulled images are cached beside the checked-out pipeline. The
`bunya_gpu` profile adds `apptainer.runOptions = '--nv'` for GPU passthrough
(dorado).
