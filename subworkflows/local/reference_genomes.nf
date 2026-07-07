/*
 * REFERENCE_GENOMES — ingest a directory of external reference genomes:
 *   normalise FASTAs -> obtain a CheckM2 quality report for them (run CheckM2, or
 *   validate a user-supplied report) -> predict proteins/nucleotides.
 *
 * Feeds three downstream consumers in the metagenome tracks:
 *   - dereplication with the HQ MAGs (COVERM_CLUSTER_HQ_REF),
 *   - read mapping to the combined set (COVERM_GENOME_HQ_REF),
 *   - the expanded gene catalogue (GENE_CATALOGUE, reference proteins).
 */

include { REFERENCE_PREP; REFERENCE_CHECKM2_VALIDATE } from '../../modules/local/reference_genomes'
include { CHECKM2_PREDICT as CHECKM2_PREDICT_REF }     from '../../modules/nf-core/checkm2/predict/main'
include { PYRODIGAL as PYRODIGAL_REFERENCE }           from '../../modules/local/pyrodigal'

workflow REFERENCE_GENOMES {
    take:
    reference_files      // collected reference genome FASTA paths
    user_checkm2         // user-supplied CheckM2 report path, or [] to run CheckM2
    checkm2_db           // CheckM2 database dir
    run_checkm2_on_refs  // bool: true => run CheckM2, false => validate user report

    main:
    ch_versions = Channel.empty()

    REFERENCE_PREP(reference_files)
    ch_genomes     = REFERENCE_PREP.out.genomes.collect()
    ch_genomes_per = REFERENCE_PREP.out.genomes.flatten().map { f -> [ [id: f.baseName], f ] }
    ch_versions    = ch_versions.mix(REFERENCE_PREP.out.versions)

    if (run_checkm2_on_refs) {
        CHECKM2_PREDICT_REF(
            ch_genomes.map { g -> [ [id: 'reference_genomes'], g ] },
            [ [:], checkm2_db ]
        )
        ch_ref_checkm2 = CHECKM2_PREDICT_REF.out.checkm2_tsv.map { m, t -> t }
        // CHECKM2_PREDICT emits versions via topic: versions (collected globally)
    } else {
        REFERENCE_CHECKM2_VALIDATE(ch_genomes, user_checkm2)
        ch_ref_checkm2 = REFERENCE_CHECKM2_VALIDATE.out.report
        ch_versions    = ch_versions.mix(REFERENCE_CHECKM2_VALIDATE.out.versions)
    }

    PYRODIGAL_REFERENCE(ch_genomes_per)
    ch_versions = ch_versions.mix(PYRODIGAL_REFERENCE.out.versions)

    emit:
    genomes     = ch_genomes                    // collected standardized FASTAs
    genomes_per = ch_genomes_per                // [ meta, fasta ] per reference
    checkm2     = ch_ref_checkm2                 // reference CheckM2 report tsv
    faa         = PYRODIGAL_REFERENCE.out.faa    // [ meta, faa ] per reference
    fna         = PYRODIGAL_REFERENCE.out.fna    // [ meta, fna ] per reference
    versions    = ch_versions
}
