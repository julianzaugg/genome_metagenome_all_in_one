# Containers

This pipeline is **containers-only** (Apptainer). Two sources:

1. **Reused nf-core modules** (`modules/nf-core/*`) ship their own biocontainer
   directive — pulled from a registry and converted by Apptainer on first use.
   Override any of them with a `withName` block in `conf/containers.config`.

2. **Bespoke local modules** use local `.sif` images under `params.container_base`.
   This defaults to `${projectDir}/containers`; profiles or the command line can
   override it for shared image locations.
   Mapped in `conf/containers.config` as `"${params.container_base}/<name>.sif"`.

### Auto-pulled from quay.io biocontainers — nothing to build
Nearly every tool is wired to a quay.io biocontainer in `conf/containers.config`
(exact build tags resolved from the quay API); Apptainer pulls and caches each on
first use. This covers — across all modes — sylph, singlem, fastplong, cleanifier,
shovill, myloasm, autocycler, polypolish, dnaapler, coverm, pyrodigal,
cd-hit, mmseqs2, dram, checkm-genome (CheckM1), nonpareil, blast (checkv
clustering), bakta, chewbbaca, parsnp, fastani, seqkit (glue), python (helpers),
and the DB downloaders. Plus the nf-core modules (fastp, spades, gtdbtk, checkm2,
geNomad, checkv).

### Must provide a local `.sif` at `params.container_base`
Tools that are not packaged on biocontainers, or whose biocontainer is not
self-contained for this pipeline, need local images:

| Image (`<name>.sif`)   | Why a local image | Used by |
|------------------------|-------------------|---------|
| `aviary_0.13.0`        | Aviary 0.13.0 requires `pixi` and prebuilt pixi environments; the quay.io biocontainer has the CLI but not `pixi` | metagenome bin recovery |
| `dorado_1.4.0`         | ONT-proprietary, not on biocontainers | Nanopore basecall/polish |
| `genomespot_1.0`       | not packaged on biocontainers | bin growth prediction (optional) |

For the **Illumina-metagenome path**, provide `aviary_0.13.0.sif` when binning is
enabled. Host removal still uses the `cleanifier` biocontainer. Supply either a
prebuilt `--cleanifier_db` index or a FASTA with `--host_ref` so the pipeline can
build the index. For Nanopore POD5 basecalling or Dorado polishing, you'll also
need `dorado`. The Dorado image used for polishing must include `samtools`,
because the wrapper sorts and indexes the Dorado aligner BAM before consensus
polishing.

> These local-`.sif` entries use `${params.container_base}`. Keep the default
> `<projectDir>/containers`, or set `--container_base` explicitly when using a
> shared image directory.

Caveats:
- **Aviary** 0.13.0 calls `pixi run` from inside its Snakemake rules. Do not use
  `quay.io/biocontainers/aviary:0.13.0--pyhdfd78af_0` for `AVIARY_RECOVER`; it
  fails with `pixi: command not found`. Build the upstream-style image, convert
  it to `aviary_0.13.0.sif`, and place it under `params.container_base`, or pass
  `--aviary_container /path/to/aviary_0.13.0.sif`.
- **CHECKV_CLUSTER** uses the same Galaxy CheckV SIF as `CHECKV_ENDTOEND` for
  blast+ plus the vendored stdlib `anicalc.py`/`aniclust.py`. The standalone
  `blast` image does not ship Python in all builds.

### Aviary SIF

The upstream Aviary Dockerfile installs `pixi`, installs `aviary-genome`, and
runs `aviary build` so the non-GPU pixi environments exist inside the image. The
important runtime property is that this command succeeds inside the SIF:

```bash
apptainer exec containers/aviary_0.13.0.sif pixi --version
apptainer exec containers/aviary_0.13.0.sif aviary --version
```

One practical build route is:

```bash
git clone https://github.com/rhysnewell/aviary /tmp/aviary
cd /tmp/aviary/docker
AVIARY_VERSION=0.13.0
sed "s/AVIARY_VERSION/$AVIARY_VERSION/g" Dockerfile.in > Dockerfile
DOCKER_BUILDKIT=1 docker build -t aviary:$AVIARY_VERSION .
apptainer build /path/to/gmaio/containers/aviary_0.13.0.sif docker-daemon://aviary:$AVIARY_VERSION
```

If your image lives elsewhere, pass `--aviary_container /path/to/aviary_0.13.0.sif`
or set that parameter in a profile.

## Apptainer config

`apptainer.enabled = true` and `autoMounts = true` are set globally. The `local`
and `bunya` profiles set `apptainer.cacheDir` to `${projectDir}/.apptainer_cache`,
so auto-pulled images are cached beside the checked-out pipeline. The
`bunya_gpu` profile adds `apptainer.runOptions = '--nv'` for GPU passthrough
(dorado).
