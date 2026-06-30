/*
 * MARKER_GENE_TREE — phylogenetic tree(s) of the user's MAGs placed alongside a
 * focused set of GTDB reference genomes, built from the GTDB-Tk marker-protein
 * alignment. Bacteria (bac120) and archaea (ar53) produce separate trees.
 *
 * Reference selection (closest-by-topology and/or related-same-order, plus an
 * optional accession list) is configured via params.marker_tree_*; with all
 * reference modes off it builds a genomes-only tree.
 */

include { MARKER_TREE_PREP        } from '../../modules/local/marker_tree'
include { MARKER_TREE_VERYFASTTREE } from '../../modules/local/marker_tree'
include { MARKER_TREE_IQTREE      } from '../../modules/local/marker_tree'

workflow MARKER_GENE_TREE {
    take:
    gtdbtk_outdir   // [ meta, gtdbtk_dir ]  (GTDBTK_CLASSIFYWF.out.gtdb_outdir)
    genomes         // collected list of user genome fastas to place
    checkm2_report  // CheckM2 quality_report.tsv (path or [])

    main:
    ch_versions = Channel.empty()

    // Build the python CLI flags from params (centralised here so the call site
    // stays clean and the module just runs the script).
    def opts = [
        "--domains bac120,ar53",
        "--min-completeness ${params.marker_tree_min_completeness}",
        "--max-contamination ${params.marker_tree_max_contamination}",
        "--closest-n ${params.marker_tree_closest_n}",
        "--related-per-order ${params.marker_tree_related_per_order}",
        params.marker_tree_use_closest ? "--use-closest" : "",
        params.marker_tree_use_related ? "--use-related" : "",
        params.marker_tree_exclude_pattern ? "--exclude-pattern '${params.marker_tree_exclude_pattern}'" : ""
    ].findAll { it }.join(' ')

    def accessions = params.marker_tree_reference_accessions
        ? file(params.marker_tree_reference_accessions, checkIfExists: true)
        : []

    MARKER_TREE_PREP(gtdbtk_outdir, checkm2_report, genomes, accessions, opts)
    ch_versions = ch_versions.mix(MARKER_TREE_PREP.out.versions)

    // One MSA per domain → [ domain, msa ]; need ≥3 taxa to build a tree.
    ch_msa = MARKER_TREE_PREP.out.msa
        .flatten()
        .map { f -> [ f.baseName.tokenize('.')[0], f ] }
        .filter { _domain, f -> f.countFasta() >= 3 }

    ch_tree = Channel.empty()
    if (params.marker_tree_builder == 'iqtree') {
        MARKER_TREE_IQTREE(ch_msa)
        ch_tree = MARKER_TREE_IQTREE.out.tree
        ch_versions = ch_versions.mix(MARKER_TREE_IQTREE.out.versions)
    } else {
        MARKER_TREE_VERYFASTTREE(ch_msa)
        ch_tree = MARKER_TREE_VERYFASTTREE.out.tree
        ch_versions = ch_versions.mix(MARKER_TREE_VERYFASTTREE.out.versions)
    }

    emit:
    tree       = ch_tree
    references = MARKER_TREE_PREP.out.references
    versions   = ch_versions
}
