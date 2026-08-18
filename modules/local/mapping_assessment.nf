/*
 * MAPPING_ASSESS — bases-mapped / percent-of-sequenced-bases metrics for a
 * CoverM read-mapping run, alongside the existing reads-mapped Count column.
 *
 * Three metrics (see bin/assess_mapping.py for the full definitions):
 *   A = full SEQ length of every primary mapped read ("how much read data
 *       mapped somewhere" -- the number to lead with for long reads)
 *   B = aligned bases of PRIMARY alignments only (CIGAR M/I/=/X)
 *   C = aligned bases of primary AND supplementary alignments -- the honest
 *       total for long reads that split across a contig/genome boundary
 *
 * Re-filters CoverM's cached BAM with `coverm filter` at the SAME identity/
 * aligned-percent thresholds as the `coverm genome` run that produced it
 * (ext.filter_args, set per-alias in conf/modules.config to match the paired
 * COVERM_GENOME / COVERM_CONTIG alias's ext.args), so these metrics describe
 * the same alignment set as the existing Count column.
 *
 * `coverm filter` evaluates its thresholds on PRIMARY records only and passes
 * every supplementary record through unconditionally (see CoverM's own
 * `coverm filter --full-help`: "Only primary, non-supplementary alignments are
 * considered"). A cleanup pass here additionally drops any supplementary
 * record whose primary did not survive the filter, so a read that fails the
 * gate contributes nothing to any metric here either -- consistent with the
 * existing Count column, at the cost of also losing that read's C-metric
 * contribution. The ONT `--min-read-aligned-percent` gate (conf/modules.config)
 * discards a chimeric/split long read outright once its primary segment alone
 * falls below the threshold; relaxing that gate would let both Count and
 * metric C credit these reads. Left as-is for now (see docs/output.md).
 */
process MAPPING_ASSESS {
    tag   { "${meta.id}:${set_label}" }
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(contig_map), path(raw_stats)
    val(set_label)

    output:
    tuple val(meta), path("${meta.id}.${set_label}.mapping_assessment.tsv"),            emit: sample_stats
    tuple val(meta), path("${meta.id}.${set_label}.mapping_assessment_per_genome.tsv"), emit: genome_stats
    path 'versions.yml', emit: versions

    script:
    def filter_args = task.ext.filter_args ?: ''
    """
    if [ -n "${filter_args}" ]; then
        coverm filter --bam-files ${bam} --output-bam-files filtered.bam \\
            ${filter_args} --threads ${task.cpus}
    else
        cp ${bam} filtered.bam
    fi

    # Drop orphan supplementary/secondary records whose primary did not survive
    # coverm filter (which only evaluates primary records against the
    # thresholds and otherwise passes every record through unchanged).
    samtools view -F 0x900 filtered.bam | cut -f1 | sort -u > keep_qnames.txt
    samtools view -h filtered.bam | awk -v OFS='\\t' '
        BEGIN { while ((getline line < "keep_qnames.txt") > 0) keep[line]=1 }
        /^@/ { print; next }
        (\$1 in keep) { print }
    ' | samtools view -@ ${task.cpus} -b -o clean.bam -

    # Strict primary-only SN numbers (A, B, denominators) -- explicit -F 0x900
    # so "bases mapped (cigar)" cannot silently include supplementary CIGAR
    # bases (samtools does this whenever the supplementary record carries a
    # real SEQ, which is minimap2's default; see assess_mapping.py docstring).
    samtools stats -@ ${task.cpus} -F 0x900 clean.bam > clean.samtools_stats.txt

    # Per-contig CIGAR pass (primary + supplementary, not secondary) for
    # metric C and the per-genome breakdown. Parses the CIGAR field directly
    # (never SEQ length), so it is correct regardless of aligner SEQ
    # conventions for supplementary records.
    samtools view -F 0x104 clean.bam | awk -F'\\t' -v OFS='\\t' '
        function cigar_bases(cig,    len, op, total) {
            total = 0
            while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
                len = substr(cig, RSTART, RLENGTH-1) + 0
                op = substr(cig, RSTART+RLENGTH-1, 1)
                if (op == "M" || op == "I" || op == "=" || op == "X") total += len
                cig = substr(cig, RSTART+RLENGTH)
            }
            return total
        }
        {
            flag = \$2; rname = \$3; cig = \$6; seqlen = length(\$10)
            is_supp = and(flag, 2048)
            bases = cigar_bases(cig)
            if (is_supp) { supp[rname]++; C[rname] += bases }
            else { prim[rname]++; A[rname] += seqlen; B[rname] += bases; C[rname] += bases }
            contigs[rname] = 1
        }
        END {
            print "contig","reads_primary","reads_supplementary","bases_A","bases_B","bases_C"
            for (r in contigs) printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n", r, prim[r]+0, supp[r]+0, A[r]+0, B[r]+0, C[r]+0
        }
    ' > per_contig.tsv

    assess_mapping.py \\
        --sample ${meta.id} --set ${set_label} \\
        --stats-file clean.samtools_stats.txt \\
        --totals-from-seqkit ${raw_stats} \\
        --per-contig per_contig.tsv \\
        --contig-map ${contig_map} \\
        --per-genome-out ${meta.id}.${set_label}.mapping_assessment_per_genome.tsv \\
        --out ${meta.id}.${set_label}.mapping_assessment.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version 2>&1 | sed 's/coverm //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    echo -e "sample\\tset\\treads_mapped\\treads_total\\tpct_reads_mapped\\treads_supplementary\\tbases_mapped_readlen\\tpct_bases_mapped_A\\tbases_mapped_cigar\\tpct_bases_mapped_B\\tbases_mapped_cigar_all\\tpct_bases_mapped_C\\tbases_total\\ttotal_source\\tnote" > ${meta.id}.${set_label}.mapping_assessment.tsv
    echo -e "${meta.id}\\t${set_label}\\t0\\t0\\t0.00\\t0\\t0\\t0.00\\t0\\t0.00\\t0\\t0.00\\t0\\tstub\\t" >> ${meta.id}.${set_label}.mapping_assessment.tsv
    echo -e "set\\tsample\\tgenome\\treads_primary\\treads_supplementary\\tbases_A\\tpct_bases_mapped_A\\tbases_B\\tpct_bases_mapped_B\\tbases_C\\tpct_bases_mapped_C" > ${meta.id}.${set_label}.mapping_assessment_per_genome.tsv
    echo '"${task.process}": {coverm: stub}' > versions.yml
    """
}
