/*
 * COMPARISON_ASSEMBLIES — pull genes from external, pre-binned assemblies
 * into the expanded gene catalogue for cross-cohort comparison. These
 * assemblies are never binned (Aviary's coverage-based recovery needs
 * matched reads, which comparison assemblies won't reliably have) -- only
 * their genes are used.
 *
 * Feeds one downstream consumer: the expanded gene catalogue (GENE_CATALOGUE,
 * alongside -- or instead of -- --reference_genomes proteins, independently
 * toggled via --reference_genomes_in_catalogue).
 */

include { COMPARISON_ASSEMBLY_PREP }                    from '../../modules/local/comparison_assemblies'
include { PYRODIGAL as PYRODIGAL_COMPARISON_ASSEMBLIES } from '../../modules/local/pyrodigal'

workflow COMPARISON_ASSEMBLIES {
    take:
    samplesheet   // path to CSV: sample, assembly

    main:
    ch_versions = Channel.empty()

    ch_rows = Channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row -> [ [id: row.sample.trim()], file(row.assembly, checkIfExists: true) ] }

    COMPARISON_ASSEMBLY_PREP(ch_rows)
    PYRODIGAL_COMPARISON_ASSEMBLIES(COMPARISON_ASSEMBLY_PREP.out.assembly)
    ch_versions = ch_versions
        .mix(COMPARISON_ASSEMBLY_PREP.out.versions)
        .mix(PYRODIGAL_COMPARISON_ASSEMBLIES.out.versions)

    emit:
    faa      = PYRODIGAL_COMPARISON_ASSEMBLIES.out.faa   // [ meta, faa ] per assembly
    fna      = PYRODIGAL_COMPARISON_ASSEMBLIES.out.fna   // [ meta, fna ] per assembly
    versions = ch_versions
}
