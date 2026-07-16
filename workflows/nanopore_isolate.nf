/*
 * NANOPORE_ISOLATE — long-read/hybrid isolate workflow.
 */

include { INPUT_CHECK }           from '../subworkflows/local/input_check'
include { LONG_READ_QC }          from '../subworkflows/local/long_read_qc'
include { ISOLATE_ANNOTATION }    from '../subworkflows/local/isolate_annotation'
include { ISOLATE_COMPARATIVE }   from '../subworkflows/local/isolate_comparative'
include { GENOME_TAXONOMY_QC }    from '../subworkflows/local/genome_taxonomy_qc'
include { MOBILE_ELEMENTS }       from '../subworkflows/local/mobile_elements'

include { FASTP }                 from '../modules/nf-core/fastp/main'
include { AUTOCYCLER_ASSEMBLE; POLYPOLISH; DNAAPLER } from '../modules/local/assembly_isolate'
include { DORADO_POLISH }         from '../modules/local/long_reads'
include { CHECKM2_PREDICT }       from '../modules/nf-core/checkm2/predict/main'
include { CHECKM1_LINEAGEWF }     from '../modules/local/checkm1'
include { COVERM_CONTIG as COVERM_CONTIG_ONT; COVERM_CONTIG as COVERM_CONTIG_SR } from '../modules/local/coverm'
include { SEQKIT_STATS }          from '../modules/local/read_stats'
include { READ_STAT_REPORT }      from '../modules/local/read_stat_report'
include { FASTQ_GZIP_TEST }       from '../modules/local/validate'
include { DUMP_SOFTWARE_VERSIONS } from '../modules/local/dump_software_versions'

def optpath = { p -> p ? file(p, checkIfExists: true) : [] }

workflow NANOPORE_ISOLATE {

    if (!params.input) { error "Mode 'nanopore_isolate' requires --input <samplesheet.csv>" }
    if (!params.skip_comparative && params.skip_annotation) {
        error "Isolate comparative analyses require Bakta GFF/FAA outputs. Disable --skip_annotation or rerun with --skip_comparative true."
    }

    INPUT_CHECK(params.input, params.mode)

    ch_direct_long = INPUT_CHECK.out.reads_long
        .filter { meta, reads -> !params.force_dorado_basecalling || !meta.has_pod5 }
    ch_pod5_basecall = INPUT_CHECK.out.pod5
        .filter { meta, pod5 -> params.force_dorado_basecalling || !meta.has_long_reads }

    LONG_READ_QC(
        ch_direct_long,
        ch_pod5_basecall,
        params.dorado_model,
        params.dorado_barcode_kit ?: '',
        params.dorado_device,
        true,
        !params.skip_porechop,
        !params.skip_qc
    )
    ch_long = LONG_READ_QC.out.reads
    ch_versions = Channel.empty()
        .mix(INPUT_CHECK.out.versions)
        .mix(LONG_READ_QC.out.versions)

    ch_short = Channel.empty()
    if (!params.skip_qc) {
        FASTP(INPUT_CHECK.out.reads_short.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
        ch_short = FASTP.out.reads
        // FASTP emits versions via topic: versions (collected globally below)
    } else {
        FASTQ_GZIP_TEST(INPUT_CHECK.out.reads_short)
        ch_short = FASTQ_GZIP_TEST.out.reads
        ch_versions = ch_versions.mix(FASTQ_GZIP_TEST.out.versions)
    }

    ch_scaffold_counts = Channel.value([])
    ch_scaffold_sr     = Channel.value([])

    ch_assembly = Channel.empty()
    if (!params.skip_assembly) {
        AUTOCYCLER_ASSEMBLE(ch_long, params.long_read_type)
        ch_assembly = AUTOCYCLER_ASSEMBLE.out.assembly
        ch_versions = ch_versions.mix(AUTOCYCLER_ASSEMBLE.out.versions)

        if (!params.skip_dorado_polish) {
            DORADO_POLISH(
                ch_assembly.join(ch_long).map { meta, assembly, reads -> [ meta, assembly, reads ] },
                params.dorado_device
            )
            ch_assembly = DORADO_POLISH.out.assembly
            ch_versions = ch_versions.mix(DORADO_POLISH.out.versions)
        }

        if (!params.skip_polypolish) {
            ch_hybrid_polish_in = ch_assembly
                .join(ch_short)
                .map { meta, assembly, reads -> [ meta, assembly, reads ] }
            POLYPOLISH(ch_hybrid_polish_in)
            ch_no_short = ch_assembly.filter { meta, assembly -> !meta.has_short_reads }
            ch_assembly = ch_no_short.mix(POLYPOLISH.out.assembly)
            ch_versions = ch_versions.mix(POLYPOLISH.out.versions)
        }

        DNAAPLER(ch_assembly)
        ch_assembly = DNAAPLER.out.assembly
        ch_versions = ch_versions.mix(DNAAPLER.out.versions)
    }

    if (!params.skip_checkm) {
        CHECKM2_PREDICT(
            ch_assembly.map { meta, fasta -> [ meta, fasta ] },
            [ [:], file(params.checkm2_db) ]
        )
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
        COVERM_CONTIG_ONT(
            ch_assembly
                .join(ch_long)
                .map { meta, scaffolds, reads -> [ meta, reads, scaffolds ] }
        )
        COVERM_CONTIG_SR(
            ch_assembly
                .join(ch_short)
                .map { meta, scaffolds, reads -> [ meta, reads, scaffolds ] }
        )
        ch_scaffold_counts = COVERM_CONTIG_ONT.out.counts.map { meta, t -> t }.collect().ifEmpty([])
        ch_scaffold_sr     = COVERM_CONTIG_SR.out.counts.map { meta, t -> t }.collect().ifEmpty([])
        ch_versions = ch_versions
            .mix(COVERM_CONTIG_ONT.out.versions)
            .mix(COVERM_CONTIG_SR.out.versions)
    }

    // --- Read-stat report (per-sample read tracking across all steps) ---
    READ_STAT_REPORT(
        'isolate',
        LONG_READ_QC.out.stats.map { meta, stage, t -> t }.collect(),
        ch_scaffold_counts,
        ch_scaffold_sr,
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
