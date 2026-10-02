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
| `dorado_1.4.0`         | ONT-proprietary, not on biocontainers (see Dorado SIF below) | Nanopore basecall/polish |
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
- **TRACS** needs 1.1.2 or newer; older builds crash on older CPUs.
  See the TRACS image section below.

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

### TRACS image

`--run_tracs` uses the `quay.io/biocontainers/tracs:1.1.4` biocontainer and needs no local image.

**Use TRACS >= 1.1.2.** Up to 1.1.1, TRACS's `setup.py` hard-coded `-march=native`, so the bioconda build carried its build host's CPU features (the 1.1.1 extension contains AVX-512).
On any older CPU it crashes:

```
Command error:
  .command.sh: line 21: 384 Illegal instruction (core dumped) tracs build-db ...
Command exit status:
  132
```

[v1.1.2](https://github.com/gtonkinhill/tracs/releases/tag/v1.1.2) dropped the flag (upstream issues #8 and #23); the target is now set only via the `TRACS_MARCH` environment variable at build time.
The 1.1.4 bioconda extension contains no AVX, BMI or SSE4.2 instructions, so it runs on any x86-64 CPU.
1.1.2 to 1.1.4 change only the build, so results match 1.1.1.

**Launch-time check.** Whatever image is used, `TRACS_PREFLIGHT` runs [`bin/tracs_smoke_test.py`](../bin/tracs_smoke_test.py), which calls every TRACS C++ kernel on tiny inputs, as soon as the pipeline starts.
An image that cannot run on the host fails within minutes, with a message pointing here, instead of days later at `TRACS_BUILD_DB`.
On a cluster it runs on one node, so it cannot vouch for nodes with different CPUs.

**Overriding.** `--tracs_container` accepts any image URI or `.sif` path, e.g. a pre-pulled copy for offline nodes:

```bash
apptainer pull containers/tracs_1.1.4.sif docker://quay.io/biocontainers/tracs:1.1.4--py312h0c6b66a_0
nextflow run . ... --run_tracs true --tracs_container containers/tracs_1.1.4.sif
```

Confirm the override took effect without launching anything:

```bash
nextflow inspect . -profile local --mode illumina_metagenome --input <samplesheet> \
    --run_tracs true --tracs_container <image> | grep -A1 TRACS_BUILD_DB
```

> Use `nextflow inspect`, not the parameter summary printed at the start of a
> run — that summary evaluates container closures against *default* parameters,
> so it shows the default image even when an override is active. The same applies
> to `--aviary_container`.

inStrain is unaffected: it is pure Python plus pysam.

### Dorado SIF

Dorado has no biocontainer and no image we can pull from Docker Hub: it is
ONT-proprietary, distributed only from ONT's own CDN under ONT's terms. Build
it locally from the checked-in recipe, which downloads the official release
tarball at build time:

```bash
apptainer build containers/dorado_1.4.0.sif containers/dorado_1.4.0.def
```

**Verify:**

```bash
apptainer exec containers/dorado_1.4.0.sif dorado --version
apptainer exec containers/dorado_1.4.0.sif samtools --version
```

The build itself needs no GPU — it only downloads and extracts the tarball —
but it does need network egress to `cdn.oxfordnanoportal.com`.

**GPU vs CPU.** The Dorado binary bundles its own CUDA runtime and only
`dlopen`s the driver's `libcuda.so` when it detects a GPU, so:
- For GPU use, the host needs an NVIDIA driver supporting CUDA driver
  ≥525.105 (for v1.4.0); GPU passthrough is already handled by
  `apptainer.runOptions = '--nv'` in `-profile bunya_gpu`.
- For CPU-only use (e.g. a local server with no GPU), set
  `--dorado_device cpu` explicitly (the default `auto` may still probe for a
  GPU). No driver, no `--nv`, and no `bunya_gpu` profile are needed in this
  case — the `local`/`bunya` profiles are sufficient.

**Rebuild for a newer version:** edit `DORADO_VERSION` in
`containers/dorado_1.4.0.def`, rebuild to a new filename
(`containers/dorado_<version>.sif`), and update the container path in
`conf/containers.config`'s `DORADO_.*` block — or point at it without editing
the pipeline via `--dorado_container /path/to/dorado_<version>.sif`.

If your image lives elsewhere, pass `--dorado_container /path/to/dorado_1.4.0.sif`
or set that parameter in a profile.

## Apptainer config

`apptainer.enabled = true` and `autoMounts = true` are set globally. The `local`
and `bunya` profiles set `apptainer.cacheDir` to `${projectDir}/.apptainer_cache`,
so auto-pulled images are cached beside the checked-out pipeline. The
`bunya_gpu` profile adds `apptainer.runOptions = '--nv'` for GPU passthrough
(dorado).
