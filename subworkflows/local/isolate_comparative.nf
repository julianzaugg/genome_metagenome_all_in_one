/*
 * ISOLATE_COMPARATIVE — build comparison groups and run group-scoped isolate
 * comparative analyses.
 */

include {
    PREP_COMPARISON_GROUPS;
    FASTANI;
    PARSNP;
    GUBBINS;
    PANAROO_RUN;
    CHEWBACCA_RUN;
    IQTREE
} from '../../modules/local/comparative'

workflow ISOLATE_COMPARATIVE {
    take:
    assemblies                 // [ meta, fasta ]
    gffs                       // [ meta, gff3 ]
    faas                       // [ meta, faa ]
    comparison_manifest        // val path string or null
    samples_include            // val path string or null
    chewbbaca_training_file    // val path string or null
    chewbbaca_cgmlst_thresholds// val
    run_fastani                // bool
    run_parsnp                 // bool
    run_panaroo                // bool
    run_chewbbaca              // bool
    run_tree                   // bool

    main:
    PREP_COMPARISON_GROUPS(
        assemblies.collect(flat: false),
        gffs.collect(flat: false),
        faas.collect(flat: false),
        comparison_manifest,
        samples_include
    )

    ch_groups = PREP_COMPARISON_GROUPS.out.groups
        .flatten()
        .map { group_dir -> [ [ id: group_dir.baseName ], group_dir ] }

    ch_fastani = Channel.empty()
    ch_parsnp = Channel.empty()
    ch_gubbins = Channel.empty()
    ch_panaroo = Channel.empty()
    ch_chewbacca = Channel.empty()
    ch_tree = Channel.empty()

    if (run_fastani) {
        FASTANI(ch_groups)
        ch_fastani = FASTANI.out.results
    }

    if (run_parsnp) {
        PARSNP(ch_groups)
        ch_parsnp = PARSNP.out.alignment
        GUBBINS(PARSNP.out.alignment)
        ch_gubbins = GUBBINS.out.polymorphic_sites
    }

    if (run_panaroo) {
        PANAROO_RUN(ch_groups)
        ch_panaroo = PANAROO_RUN.out.core_alignment
        if (run_tree) {
            IQTREE(PANAROO_RUN.out.core_alignment)
            ch_tree = IQTREE.out.tree
        }
    }

    if (run_chewbbaca) {
        CHEWBACCA_RUN(ch_groups, chewbbaca_training_file, chewbbaca_cgmlst_thresholds)
        ch_chewbacca = CHEWBACCA_RUN.out.outdir
    }

    emit:
    groups    = PREP_COMPARISON_GROUPS.out.groups
    fastani   = ch_fastani
    parsnp    = ch_parsnp
    gubbins   = ch_gubbins
    panaroo   = ch_panaroo
    chewbacca = ch_chewbacca
    tree      = ch_tree
}
