# Containers

This pipeline is **containers-only** (Apptainer). Two sources:

1. **Reused nf-core modules** (`modules/nf-core/*`) ship their own biocontainer
   directive — pulled from a registry and converted by Apptainer on first use.
   Override any of them with a `withName` block in `conf/containers.config`.

2. **Bespoke local modules** use local `.sif` images under `params.container_base`
   (set per profile: `/scratch/project/a_ace/containers` on Bunya,
   `/srv/db/containers` locally). Mapped in `conf/containers.config` as
   `"${params.container_base}/<name>.sif"`.

## .sif image checklist (bespoke tools)

Provide/confirm these images at `params.container_base` (filenames are editable
in `conf/containers.config`):

| Image (`<name>.sif`)   | Tool(s)                          | Used by mode |
|------------------------|----------------------------------|--------------|
| `sylph_0.9.0`          | sylph                            | metagenomes |
| `singlem_0.19.0`       | singlem                          | metagenomes |
| `cleanifier_1.3.0`     | cleanifier                       | host removal |
| `minimap2_2.28`        | minimap2 + samtools              | host removal / mapping |
| `aviary_0.12.0`        | Aviary (complex deps)            | metagenomes |
| `coverm_0.7.0`         | CoverM                           | derep + mapping |
| `pyrodigal_3.6.3`      | pyrodigal                        | gene calling |
| `cdhit_4.8.1`          | cd-hit                           | gene catalogue |
| `seqkit_2.8`           | small fasta munging              | prep/process steps |
| `python_3.11`          | bin/ helper scripts              | tabulate, etc. |
| `dram_1.4.6`           | DRAM                             | annotation |
| `checkm1_1.2.3`        | CheckM1                          | bin QC |
| `nonpareil_3.4.1`      | nonpareil                        | coverage |
| `genomespot_1.0`       | GenomeSPOT                       | bin growth prediction |
| `checkv_1.0.3`         | CheckV + blast (anicalc/aniclust)| clustering |
| `fastplong_0.4.1`      | fastplong                        | nanopore QC |
| `myloasm_0.5.1`        | myloasm                          | nanopore meta assembly |
| `autocycler_0.6.1`     | autocycler                       | nanopore isolate assembly |
| `dorado_1.4.0`         | dorado (GPU)                     | basecall/polish |
| `polypolish_0.6.1`     | polypolish                       | hybrid polish |
| `dnaapler_1.3.0`       | dnaapler                         | reorientation |
| `shovill_1.1.0`        | shovill                          | illumina isolate assembly |
| `chewbbaca_3.3.10`     | chewBBACA                        | cgMLST |
| `parsnp_2.1.4`         | parsnp                           | core alignment |
| `fastani_1.33`         | fastANI                          | comparative |

## Apptainer config

`apptainer.enabled = true` and `autoMounts = true` are set globally;
`apptainer.cacheDir` is set per profile. The `bunya_gpu` profile adds
`apptainer.runOptions = '--nv'` for GPU passthrough (dorado).
