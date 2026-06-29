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

The image is maintained at
[julianzaugg/aviary](https://github.com/julianzaugg/aviary) and built
automatically by GitHub Actions on every `v*` tag push or manual trigger.
It installs pixi, clones the upstream aviary repo, pip-installs
`aviary-genome`, runs `aviary build` to pre-build all non-GPU pixi
environments, and bakes the conda env PATH into a static `/entrypoint.sh`
so aviary starts without pixi at the outer level.

**Pull from Docker Hub (recommended):**

```bash
# On the machine where the SIF will be used (e.g. page, Bunya)
apptainer pull docker://julianzaugg/aviary:0.13.0
mv aviary_0.13.0.sif /path/to/gmaio/containers/
```

**Verify:**

```bash
apptainer run containers/aviary_0.13.0.sif --help
```

**Rebuild the Docker image** (e.g. to update the aviary version):

1. Update the `git checkout` line and version in
   [julianzaugg/aviary/docker/Dockerfile](https://github.com/julianzaugg/aviary/blob/main/docker/Dockerfile)
2. Commit and push to `main`
3. Go to Actions → "Build and push Docker image" → Run workflow → enter the new tag
4. Re-pull the SIF on the compute host

**Runtime notes:**

- Aviary's Snakemake rules call `pixi run` internally for subtools (coverm,
  rosella, etc.). On servers with NFS home directories, set `PIXI_CACHE_DIR`
  to local scratch — the module already exports
  `PIXI_CACHE_DIR=/tmp/pixi-cache-$USER` automatically.
- If pixi still hits read-only filesystem errors during a real run, add
  `runOptions = '--writable-tmpfs'` to the `apptainer {}` block in your
  local config.
- The SIF requires `procps` (`ps`) for Nextflow task metrics. Images built
  from commit `21db7ff` onward include it. For older SIFs, add
  `runOptions = '--bind /usr/bin/ps:/usr/bin/ps'` to the `apptainer {}` block.

If your image lives elsewhere, pass `--aviary_container /path/to/aviary_0.13.0.sif`
or set that parameter in a profile.

## Apptainer config

`apptainer.enabled = true` and `autoMounts = true` are set globally. The `local`
and `bunya` profiles set `apptainer.cacheDir` to `${projectDir}/.apptainer_cache`,
so auto-pulled images are cached beside the checked-out pipeline. The
`bunya_gpu` profile adds `apptainer.runOptions = '--nv'` for GPU passthrough
(dorado).
