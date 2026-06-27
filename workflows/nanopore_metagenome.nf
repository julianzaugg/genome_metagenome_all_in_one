/*
 * NANOPORE_METAGENOME — long-read metagenome workflow.
 */

include { INPUT_CHECK }        from '../subworkflows/local/input_check'
include { LONG_READ_QC }       from '../subworkflows/local/long_read_qc'
include { READ_PROFILING }     from '../subworkflows/local/read_profiling'
include { HOST_REMOVAL }       from '../subworkflows/local/host_removal'
include { GENE_CATALOGUE }     from '../subworkflows/local/gene_catalogue'
include { GENOME_TAXONOMY_QC } from '../subworkflows/local/genome_taxonomy_qc'
include { MOBILE_ELEMENTS }    from '../subworkflows/local/mobile_elements'

include { MYLOASM }                     from '../modules/local/assembly_isolate'
include { DORADO_POLISH }               from '../modules/local/long_reads'
include { CHECKM2_PREDICT }             from '../modules/nf-core/checkm2/predict/main'
include { AVIARY_RECOVER; AVIARY_COLLECT_BINS } from '../modules/local/aviary'
include { COVERM_CLUSTER; COVERM_GENOME as COVERM_GENOME_ONT; COVERM_CONTIG as COVERM_CONTIG_ONT } from '../modules/local/coverm'
include { CHECKM1_LINEAGEWF }           from '../modules/local/checkm1'
include { PYRODIGAL as PYRODIGAL_SCAFFOLDS } from '../modules/local/pyrodigal'
include { NONPAREIL }                   from '../modules/local/nonpareil'
include { CLEANIFIER_INDEX }            from '../modules/local/host_removal'

def optpath = { p -> p ? file(p, checkIfExists: true) : [] }

workflow NANOPORE_METAGENOME {

    if (!params.input) { error "Mode 'nanopore_metagenome' requires --input <samplesheet.csv>" }

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
    ch_reads = LONG_READ_QC.out.reads

    ch_sylph_tax_meta = params.sylph_tax_metadata
        ? Channel.fromPath(params.sylph_tax_metadata, checkIfExists: true).collect()
        : Channel.value([])
    READ_PROFILING(
        ch_reads,
        Channel.fromPath(params.sylph_db, checkIfExists: true).collect(),
        ch_sylph_tax_meta,
        file(params.singlem_metapackage, checkIfExists: true),
        !params.skip_sylph,
        !params.skip_sylph && params.sylph_tax_metadata != null,
        !params.skip_singlem
    )

    if (!params.skip_host_removal) {
        if (params.cleanifier_db) {
            ch_cleanifier_index = Channel.value(file(params.cleanifier_db, checkIfExists: true))
        } else {
            if (!params.host_ref) {
                error "Host removal is enabled but neither --cleanifier_db nor --host_ref is set. Provide a Cleanifier .filter index, provide a FASTA with --host_ref to build one, or rerun with --skip_host_removal true."
            }
            CLEANIFIER_INDEX(file(params.host_ref, checkIfExists: true), params.cleanifier_nobjects ?: '')
            ch_cleanifier_index = CLEANIFIER_INDEX.out.index
        }
        HOST_REMOVAL(ch_reads, ch_cleanifier_index)
        ch_clean = HOST_REMOVAL.out.reads
    } else {
        ch_clean = ch_reads
    }

    ch_assembly = Channel.empty()
    if (!params.skip_assembly) {
        MYLOASM(ch_clean)
        ch_assembly = MYLOASM.out.assembly

        if (!params.skip_dorado_polish) {
            DORADO_POLISH(
                ch_assembly.join(ch_clean).map { meta, assembly, reads -> [ meta, assembly, reads ] },
                params.dorado_device
            )
            ch_assembly = DORADO_POLISH.out.assembly
        }
    }

    if (!params.skip_assembly && !params.skip_read_mapping) {
        COVERM_CONTIG_ONT(
            ch_assembly
                .join(ch_clean)
                .map { meta, scaffolds, reads -> [ meta, reads, scaffolds ] }
        )
    }

    if (!params.skip_binning) {
        AVIARY_RECOVER(
            ch_assembly.join(ch_clean),
            file(params.gtdbtk_db,  checkIfExists: true),
            file(params.checkm2_db, checkIfExists: true),
            file(params.eggnog_db,  checkIfExists: true)
        )
        AVIARY_COLLECT_BINS(AVIARY_RECOVER.out.bins.map { meta, bins -> bins }.flatten().collect())
        ch_all_bins = AVIARY_COLLECT_BINS.out.bins

        ch_checkm2_tsv = Channel.value([])
        ch_checkm1_tsv = Channel.value([])
        if (!params.skip_checkm) {
            CHECKM2_PREDICT(
                ch_all_bins.map { bins -> [ [id: 'all_bins'], bins ] },
                [ [:], file(params.checkm2_db) ]
            )
            ch_checkm2_tsv = CHECKM2_PREDICT.out.checkm2_tsv.map { m, t -> t }
        }
        if (params.run_checkm1) {
            CHECKM1_LINEAGEWF(ch_all_bins)
            ch_checkm1_tsv = CHECKM1_LINEAGEWF.out.summary
        }

        if (!params.skip_dereplication) {
            COVERM_CLUSTER(ch_all_bins, ch_checkm2_tsv, ch_checkm1_tsv)
            ch_reps      = COVERM_CLUSTER.out.representatives.collect()
            ch_per_rep   = COVERM_CLUSTER.out.representatives.flatten()
                                .map { f -> [ [id: f.baseName], f ] }
        } else {
            ch_reps    = ch_all_bins
            ch_per_rep = AVIARY_COLLECT_BINS.out.bins.flatten().map { b -> [ [id: b.baseName], b ] }
        }

        if (!params.skip_read_mapping) {
            COVERM_GENOME_ONT(ch_clean, ch_reps)
        }

        GENOME_TAXONOMY_QC(
            ch_reps,
            ch_per_rep,
            file(params.gtdbtk_db, checkIfExists: true),
            optpath(params.genomespot_models),
            !params.skip_taxonomy,
            params.run_genomespot,
            params.run_barrnap
        )
    }

    if (!params.skip_gene_catalogue) {
        PYRODIGAL_SCAFFOLDS(ch_assembly)
        GENE_CATALOGUE(
            PYRODIGAL_SCAFFOLDS.out.faa,
            PYRODIGAL_SCAFFOLDS.out.fna,
            params.catalogue_identities,
            optpath(params.dram_db),
            !params.skip_annotation
        )
    }

    if (!params.skip_mobile_elements) {
        MOBILE_ELEMENTS(
            ch_assembly,
            file(params.genomad_db, checkIfExists: true),
            file(params.checkv_db,  checkIfExists: true)
        )
    }

    if (params.run_nonpareil) {
        NONPAREIL(ch_clean)
    }
}
