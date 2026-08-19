/*
 * STRAIN_COMPARISON — do different samples carry the SAME strain, or different strains
 * of the same species?
 *
 * Dereplication at 95% ANI only establishes that two samples share a *species*. To go
 * below that, every sample's reads are mapped to ONE shared cross-sample reference set
 * (the dereplicated MAGs), and the resulting per-sample allele profiles are compared
 * pairwise. A shared reference is essential: per-sample references would give each
 * sample its own coordinate system and nothing would be comparable.
 *
 * Two independent tools, because neither covers both platforms well:
 *   inStrain — the established popANI approach. Short reads only: its SNV model assumes
 *              low per-read error, and no published parameterisation makes it sound for
 *              ONT (swapping in a long-read aligner does not change that).
 *   TRACS    — empirical Bayes over variable coverage, emits pairwise SNP distances and
 *              transmission clusters, and supports both short and long reads.
 *
 * The 95% ANI clustering upstream is load-bearing here: one reference per species cluster
 * means every sample's reads for a species land on the SAME reference. Clustering tighter
 * would split strains of one species across references and they would never be compared.
 */

include { STRAIN_GENOME_FILTER  } from '../../modules/local/strain_comparison'
include { STRAIN_REFERENCE_PREP } from '../../modules/local/strain_comparison'
include { BOWTIE2_STRAIN_BUILD  } from '../../modules/local/strain_comparison'
include { BOWTIE2_STRAIN_ALIGN  } from '../../modules/local/strain_comparison'
include { INSTRAIN_PROFILE      } from '../../modules/local/strain_comparison'
include { INSTRAIN_COMPARE      } from '../../modules/local/strain_comparison'
include { INSTRAIN_SUMMARISE    } from '../../modules/local/strain_comparison'
include { TRACS_BUILD_DB        } from '../../modules/local/strain_comparison'
include { TRACS_ALIGN           } from '../../modules/local/strain_comparison'
include { TRACS_ALIGN as TRACS_ALIGN_ONT } from '../../modules/local/strain_comparison'
include { TRACS_COMBINE         } from '../../modules/local/strain_comparison'
include { TRACS_DISTANCE        } from '../../modules/local/strain_comparison'
include { TRACS_CLUSTER         } from '../../modules/local/strain_comparison'

workflow STRAIN_COMPARISON {
    take:
    reads           // [ meta, reads ]  — QC'd / host-removed reads, one entry per sample
    genomes         // collected list of representative genome fastas (strain_genome_source set)
    checkm2_report  // CheckM2 quality_report.tsv (path or [])
    checkm1_report  // CheckM1 summary tsv (path or [])
    run_instrain    // bool
    run_tracs       // bool
    long_reads      // bool — selects the TRACS minimap2 preset

    main:
    ch_versions = Channel.empty()

    // Restrict to genomes good enough for strain calling. Both tools consume the SAME
    // filtered set, so their calls stay directly comparable.
    STRAIN_GENOME_FILTER(genomes, checkm2_report, checkm1_report)
    ch_versions = ch_versions.mix(STRAIN_GENOME_FILTER.out.versions)

    // The filter can legitimately empty the set (e.g. no MAG reaches 90/5). Drop the
    // empty emission rather than letting downstream tools fail on an unmatched glob --
    // the same tolerance COVERM_GENOME applies to an empty HQ set.
    ch_ref_genomes = STRAIN_GENOME_FILTER.out.genomes
        .collect()
        .filter { gs -> gs && gs.size() > 0 }

    ch_instrain_compare = Channel.empty()
    ch_instrain_summary = Channel.empty()
    ch_tracs_distances  = Channel.empty()
    ch_tracs_clusters   = Channel.empty()

    if (run_instrain) {
        // Combined FASTA + .stb, contigs prefixed with their bin name (bins from separate
        // per-sample assemblies can otherwise share contig names).
        STRAIN_REFERENCE_PREP(ch_ref_genomes)
        BOWTIE2_STRAIN_BUILD(STRAIN_REFERENCE_PREP.out.fasta)
        BOWTIE2_STRAIN_ALIGN(reads, BOWTIE2_STRAIN_BUILD.out.index)

        INSTRAIN_PROFILE(
            BOWTIE2_STRAIN_ALIGN.out.bam,
            STRAIN_REFERENCE_PREP.out.fasta,
            STRAIN_REFERENCE_PREP.out.stb
        )
        INSTRAIN_COMPARE(
            INSTRAIN_PROFILE.out.profile.map { _meta, p -> p }.collect(),
            STRAIN_REFERENCE_PREP.out.stb
        )
        INSTRAIN_SUMMARISE(INSTRAIN_COMPARE.out.compare)

        ch_instrain_compare = INSTRAIN_COMPARE.out.compare
        ch_instrain_summary = INSTRAIN_SUMMARISE.out.summary
        ch_versions = ch_versions
            .mix(STRAIN_REFERENCE_PREP.out.versions)
            .mix(BOWTIE2_STRAIN_BUILD.out.versions)
            .mix(BOWTIE2_STRAIN_ALIGN.out.versions)
            .mix(INSTRAIN_PROFILE.out.versions)
            .mix(INSTRAIN_COMPARE.out.versions)
            .mix(INSTRAIN_SUMMARISE.out.versions)
    }

    if (run_tracs) {
        // A build-db database embeds both the sourmash index and the genomes themselves,
        // so `tracs align` never falls back to downloading GTDB references from Genbank.
        TRACS_BUILD_DB(ch_ref_genomes)
        ch_versions = ch_versions.mix(TRACS_BUILD_DB.out.versions)

        if (long_reads) {
            TRACS_ALIGN_ONT(reads, TRACS_BUILD_DB.out.db)
            ch_tracs_align = TRACS_ALIGN_ONT.out.alignment
            ch_versions = ch_versions.mix(TRACS_ALIGN_ONT.out.versions)
        } else {
            TRACS_ALIGN(reads, TRACS_BUILD_DB.out.db)
            ch_tracs_align = TRACS_ALIGN.out.alignment
            ch_versions = ch_versions.mix(TRACS_ALIGN.out.versions)
        }

        // align writes one posterior-count fasta per (sample, reference); combine
        // transposes that into one MSA per reference across samples, which is the only
        // form `distance` accepts.
        TRACS_COMBINE(ch_tracs_align.map { _meta, d -> d }.collect())
        TRACS_DISTANCE(TRACS_COMBINE.out.combined)
        TRACS_CLUSTER(TRACS_DISTANCE.out.distances)

        ch_tracs_distances = TRACS_DISTANCE.out.distances
        ch_tracs_clusters  = TRACS_CLUSTER.out.clusters
        ch_versions = ch_versions
            .mix(TRACS_COMBINE.out.versions)
            .mix(TRACS_DISTANCE.out.versions)
            .mix(TRACS_CLUSTER.out.versions)
    }

    emit:
    reference_genomes = ch_ref_genomes
    reference_report  = STRAIN_GENOME_FILTER.out.report
    instrain_compare  = ch_instrain_compare
    instrain_summary  = ch_instrain_summary
    tracs_distances   = ch_tracs_distances
    tracs_clusters    = ch_tracs_clusters
    versions          = ch_versions
}
