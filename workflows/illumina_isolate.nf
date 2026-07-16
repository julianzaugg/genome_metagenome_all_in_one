/*
 * ILLUMINA_ISOLATE — short-read isolate workflow.
 */

include { INPUT_CHECK }           from '../subworkflows/local/input_check'
include { ISOLATE_ANNOTATION }    from '../subworkflows/local/isolate_annotation'
include { ISOLATE_COMPARATIVE }   from '../subworkflows/local/isolate_comparative'
include { GENOME_TAXONOMY_QC }    from '../subworkflows/local/genome_taxonomy_qc'
include { MOBILE_ELEMENTS }       from '../subworkflows/local/mobile_elements'

include { FASTP }                  from '../modules/nf-core/fastp/main'
include { SHOVILL }                from '../modules/local/assembly_isolate'
include { CHECKM2_PREDICT }        from '../modules/nf-core/checkm2/predict/main'
include { CHECKM1_LINEAGEWF }      from '../modules/local/checkm1'
include { COVERM_CONTIG }          from '../modules/local/coverm'
include { SEQKIT_STATS }           from '../modules/local/read_stats'
include { READ_STAT_REPORT }       from '../modules/local/read_stat_report'
include { FASTQ_GZIP_TEST }        from '../modules/local/validate'
include { DUMP_SOFTWARE_VERSIONS } from '../modules/local/dump_software_versions'

def optpath = { p -> p ? file(p, checkIfExists: true) : [] }

workflow ILLUMINA_ISOLATE {

    if (!params.input) { error "Mode 'illumina_isolate' requires --input <samplesheet.csv>" }
    if (!params.skip_assembly && params.skip_qc) {
        log.warn "[gmaio] Illumina isolate assembly will use raw reads because --skip_qc true."
    }
    if (!params.skip_comparative && params.skip_annotation) {
        error "Isolate comparative analyses require Bakta GFF/FAA outputs. Disable --skip_annotation or rerun with --skip_comparative true."
    }

    INPUT_CHECK(params.input, params.mode)
    FASTQ_GZIP_TEST(INPUT_CHECK.out.reads_short)
    ch_reads = FASTQ_GZIP_TEST.out.reads
    ch_read_stats = ch_reads.map { meta, reads -> [ meta, 'raw', reads ] }
    ch_versions = Channel.empty()
        .mix(INPUT_CHECK.out.versions)
        .mix(FASTQ_GZIP_TEST.out.versions)

    if (!params.skip_qc) {
        FASTP(ch_reads.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
        ch_qc = FASTP.out.reads
        ch_read_stats = ch_read_stats.mix(FASTP.out.reads.map { meta, reads -> [ meta, 'fastp', reads ] })
        // FASTP emits versions via topic: versions (collected globally below)
    } else {
        ch_qc = ch_reads
    }
    SEQKIT_STATS(ch_read_stats)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions)

    ch_scaffold_counts = Channel.value([])

    ch_assembly = Channel.empty()
    if (!params.skip_assembly) {
        SHOVILL(ch_qc)
        ch_assembly = SHOVILL.out.assembly
        ch_versions = ch_versions.mix(SHOVILL.out.versions)
    }

    ch_checkm2_tsv = Channel.empty()
    if (!params.skip_checkm) {
        CHECKM2_PREDICT(
            ch_assembly.map { meta, fasta -> [ meta, fasta ] },
            [ [:], file(params.checkm2_db) ]
        )
        ch_checkm2_tsv = CHECKM2_PREDICT.out.checkm2_tsv
        // CHECKM2_PREDICT emits versions via topic: versions (collected globally below)
    }
    if (params.run_checkm1) {
        CHECKM1_LINEAGEWF(ch_assembly.map { meta, fasta -> fasta }.collect())
        ch_versions = ch_versions.mix(CHECKM1_LINEAGEWF.out.versions)
    }

    GENOME_TAXONOMY_QC(
        ch_assembly.map { meta, fasta -> fasta }.collect(),
        ch_assembly,
        file(params.gtdbtk_db, checkIfExists: true),
        optpath(params.genomespot_models),
        optpath(params.dram_db),
        !params.skip_taxonomy,
        params.run_genomespot,
        params.run_barrnap,
        params.run_dram_bins,
        Channel.value([]),   // no external reference genomes in isolate mode
        false                // no USERREF_ prefix to restore
    )
    ch_versions = ch_versions.mix(GENOME_TAXONOMY_QC.out.versions)

    ISOLATE_ANNOTATION(
        ch_assembly,
        file(params.bakta_db, checkIfExists: true),
        optpath(params.bakta_reference_proteins),
        optpath(params.amrfinder_db),
        !params.skip_annotation,
        !params.skip_annotation,
        !params.skip_annotation
    )
    ch_versions = ch_versions.mix(ISOLATE_ANNOTATION.out.versions)

    if (!params.skip_mobile_elements) {
        MOBILE_ELEMENTS(
            ch_assembly,
            file(params.genomad_db, checkIfExists: true),
            file(params.checkv_db,  checkIfExists: true)
        )
        ch_versions = ch_versions.mix(MOBILE_ELEMENTS.out.versions)
    }

    if (!params.skip_read_mapping) {
        COVERM_CONTIG(
            ch_assembly
                .join(ch_qc)
                .map { meta, scaffolds, reads -> [ meta, reads, scaffolds ] }
        )
        ch_scaffold_counts = COVERM_CONTIG.out.counts.map { meta, t -> t }.collect().ifEmpty([])
        ch_versions = ch_versions.mix(COVERM_CONTIG.out.versions)
    }

    // --- Read-stat report (per-sample read tracking across all steps) ---
    READ_STAT_REPORT(
        'isolate',
        SEQKIT_STATS.out.stats.map { meta, stage, t -> t }.collect(),
        ch_scaffold_counts,
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(READ_STAT_REPORT.out.versions)

    if (!params.skip_comparative) {
        ISOLATE_COMPARATIVE(
            ch_assembly,
            ISOLATE_ANNOTATION.out.gff,
            ISOLATE_ANNOTATION.out.faa,
            params.comparison_manifest ?: '',
            params.samples_include ?: '',
            params.chewbbaca_training_file ?: '',
            params.chewbbaca_cgmlst_thresholds,
            params.run_fastani,
            params.run_parsnp,
            params.run_panaroo,
            params.run_chewbbaca,
            params.run_tree
        )
        ch_versions = ch_versions.mix(ISOLATE_COMPARATIVE.out.versions)
    }

    // --- Software versions manifest ---
    ch_nfcore_versions = channel.topic('versions')
        .collectFile(name: 'nfcore_versions.yml', newLine: true) { process, tool, version ->
            "\"${process}\":\n    ${tool}: ${version}"
        }
    DUMP_SOFTWARE_VERSIONS(ch_versions.mix(ch_nfcore_versions).collect())
}
