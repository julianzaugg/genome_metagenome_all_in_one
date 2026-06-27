# Changelog

All notable changes to this pipeline are documented here.
Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- Initial scaffold: four `--mode` entry workflows, nf-core-style layout.
- `illumina_metagenome` workflow built end-to-end (fastp → sylph + singlem →
  cleanifier → metaspades → Aviary → CoverM dereplication/mapping → pyrodigal →
  cd-hit gene catalogue → DRAM → GTDB-Tk → CheckM1/CheckM2 → nonpareil →
  genomespot → geNomad → CheckV → ANI clustering).
- `nanopore_metagenome`, `illumina_isolate`, `nanopore_isolate` implemented as
  production DSL2 workflows with stub-test coverage.
- Nanopore long-read QC/basecalling subworkflow with optional Dorado basecalling
  for POD5-only rows or forced re-basecalling.
- Isolate annotation and group-scoped comparative subworkflows, including mixed
  sample/reference comparison manifests and overlapping comparison groups.
- Containers-only config; profiles `bunya`, `bunya_gpu`, `local`, `test`.
- Reference-database params + optional `download_dbs` entry workflow.
