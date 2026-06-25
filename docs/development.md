# Development

## Running locally for wiring checks

The pipeline is containers-only, so a *real* run needs Apptainer (HPC/Linux).
But you can validate **wiring/topology** anywhere with `-stub` + the self-contained
`test` profile (no tools, no containers, tiny resources):

```bash
nextflow run . -profile test --mode illumina_metagenome -stub
nextflow run . -profile test --mode download_dbs --download_db checkm2,genomad -stub
```

`-profile test` sets the local executor, disables Apptainer, points all `--*_db`
params at `tests/data/dbs/`, and puts `tests/shims/` on `PATH` (tiny fake
`fastp`/`spades.py`/`gtdbtk`/… that satisfy the nf-core modules' `eval()`
version-capture, which runs even under `-stub`).

## Nextflow / Java on macOS (this repo's dev box)

The system `nextflow` (23.10.1, capsule launcher) fails to parse the installed
Zulu JDK 21 version string. Workaround used here (not committed): the official
non-capsule launcher pinned to 24.10.5 with `JAVA_HOME` set to Zulu 21. The
`nf-schema` plugin must be present in `~/.nextflow/plugins/nf-schema-2.1.1/`.
On Bunya/the Linux server, use the site Nextflow + Apptainer as normal.

## Managing nf-core modules

```bash
nf-core modules list local
nf-core modules install <tool>
nf-core modules update <tool>
```
`nf-core` shells out to `nextflow config`, so the working `nextflow` must be on
`PATH`. Reused-module containers are overridden in `conf/containers.config`.

## Tests

`nf-test` definitions live in `tests/`. Run with `nf-test test` (requires
`nf-test` installed). The stub run above is the lightweight integration check.

## Adding a step

1. Add/installed the module (`modules/local/*` or `nf-core modules install`).
2. Map its container + `ext.args`/`publishDir` in `conf/containers.config` and
   `conf/modules.config`.
3. Wire it into the relevant `subworkflows/local/*` and `workflows/<mode>.nf`,
   guarded by a `--skip_*`/`--run_*` param.
4. Give it a `stub:` block and re-run the stub check.
