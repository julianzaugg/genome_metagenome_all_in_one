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
| `tracs_1.1.1`          | **only if the biocontainer SIGILLs on your CPU** — upstream compiles with `-march=native`, so the published image is not portable (see below) | strain comparison (`--run_tracs`) |

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
- **TRACS** works from the quay.io biocontainer on many hosts, but not all — see
  the TRACS SIF section below. It is wired to the biocontainer by default and
  only needs a local image if it crashes with exit status 132.

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

### TRACS SIF (only if the biocontainer crashes)

`--run_tracs` uses `quay.io/biocontainers/tracs` by default and needs no local
image. But if `TRACS_BUILD_DB` (or any `TRACS_*` process) dies like this:

```
Command error:
  .command.sh: line 2: 51 Illegal instruction (core dumped) tracs build-db ...
Command exit status:
  132
```

that is `SIGILL` — the binary uses CPU instructions your host does not have.

**Cause.** TRACS's `setup.py` hard-codes

```python
extra_compile_args = ["-O3", "-ffast-math", "-march=native"]
```

so its pybind11 extension is compiled for whatever CPU built it. bioconda's
build host is therefore baked into the published image: it runs on CPUs at
least as new as that machine, and SIGILLs on anything older. Upstream issue
[#8](https://github.com/gtonkinhill/tracs/issues/8) is exactly this failure and
is *closed*, but no fix ever landed — the reporter closed it after a bioconda
rebuild happened to suit their hardware. The flag is still present on `main`
(v1.1.1, the newest release and newest build). So there is no version to
upgrade to, trying other build tags is a coin flip, and an image that works
today can break on a future rebuild.

**Fix — build it on the machine that will run it**, so `-march=native` targets
your own CPU:

```bash
# On the compute host (page, Bunya, ...) -- NOT on a newer machine
apptainer build containers/tracs_1.1.1.sif containers/tracs_1.1.1.def
```

The build takes a few minutes (~380 MB of conda packages plus the C++
extension). It verifies `import TRACS`, `tracs --version`, every subcommand's
`--help`, and the presence of samtools/minimap2/htsbox/sourmash before
finishing — so a bad image fails at build time rather than hours into a
pipeline run. The `import` is the check that matters: that is what SIGILLs when
the extension does not match the CPU.

`cxx-compiler` is in the package list deliberately. Without a C++ toolchain the
build fails late inside pybind11 with a misleading message:

```
RuntimeError: Unsupported compiler -- at least C++11 support is needed!
```

conda-forge's compilers also ship activation scripts that do **not** run inside
`%post`, so the definition file exports `CC`/`CXX` explicitly.

Then point the pipeline at it:

```bash
nextflow run . ... --run_tracs true --tracs_container /path/to/tracs_1.1.1.sif
```

Confirm it took effect without launching anything:

```bash
nextflow inspect . -profile local --mode illumina_metagenome --input <samplesheet> \
    --run_tracs true --tracs_container /path/to/tracs_1.1.1.sif | grep -A1 TRACS_BUILD_DB
```

> Use `nextflow inspect`, not the parameter summary printed at the start of a
> run — that summary evaluates container closures against *default* parameters,
> so it shows the biocontainer even when an override is active. The same applies
> to `--aviary_container`.

**Mixed hardware.** A SIF built on a newer CPU reintroduces the same crash on
older nodes. For a heterogeneous cluster, build once against a portable
baseline instead of `native` by editing the marked line near the top of
`%post` in `containers/tracs_1.1.1.def`:

```bash
TRACS_MARCH="x86-64-v2"    # SSE4.2 baseline (~2009+); x86-64-v3 for AVX2 (~2013+)
```

It is a plain shell variable rather than a build argument so it works on any
Apptainer version. The `.def` applies it by patching `setup.py`, because
setuptools appends `extra_compile_args` last and `CFLAGS` cannot override them.

**Diagnosing.** To confirm SIGILL is the extension rather than a dependency:

```bash
apptainer exec <image> python -c "import TRACS; print('extension loads OK')"
lscpu | grep -oE 'avx[0-9a-z_]+' | sort -u     # what your CPU actually supports
```

inStrain is unaffected — it is pure Python plus pysam. If TRACS blocks you, run
with `--run_instrain true --run_tracs false` and add TRACS once the image is
built.

## Apptainer config

`apptainer.enabled = true` and `autoMounts = true` are set globally. The `local`
and `bunya` profiles set `apptainer.cacheDir` to `${projectDir}/.apptainer_cache`,
so auto-pulled images are cached beside the checked-out pipeline. The
`bunya_gpu` profile adds `apptainer.runOptions = '--nv'` for GPU passthrough
(dorado).
