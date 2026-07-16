/*
 * External reference genomes — normalise a user-supplied directory of genome
 * FASTAs and (optionally) validate a user-supplied CheckM2 report for them.
 *
 * The standardized files keep the ORIGINAL basename stem (only the extension is
 * normalised to .fasta and .gz inputs are decompressed) so that a user-supplied
 * CheckM2 quality_report.tsv — whose "Name" column holds the original stems —
 * still matches. Collision-safety against MAG names is enforced later, at
 * COVERM_CLUSTER_HQ_REF, where refs and bins are staged into one folder.
 */

process REFERENCE_PREP {
    label 'process_low'

    input:
    path(refs, stageAs: 'refs_in/*')

    output:
    path 'reference_genomes/*.fasta', emit: genomes
    path 'gtdbtk_refs/*.fasta',       emit: gtdbtk_genomes
    path 'versions.yml',              emit: versions

    script:
    // gtdbtk_refs/ holds a copy of each reference whose stem is prefixed with the
    // reserved token USERREF_ so GTDB-Tk sees a non-accession id — avoids the hard
    // error GTDB-Tk raises when a user genome id (e.g. GCA_/GCF_...) collides with a
    // GTDB reference accession. The prefix is stripped again in the tree labels and
    // the published GTDB-Tk summary.
    """
    mkdir -p reference_genomes gtdbtk_refs
    : > seen_stems.txt
    for f in refs_in/*; do
        [ -e "\$f" ] || continue
        b=\$(basename "\$f")
        # strip a trailing .gz, then a known FASTA extension, to get the stem
        stem="\$b"
        case "\$stem" in *.gz) stem="\${stem%.gz}" ;; esac
        case "\$stem" in
            *.fasta) stem="\${stem%.fasta}" ;;
            *.fa)    stem="\${stem%.fa}" ;;
            *.fna)   stem="\${stem%.fna}" ;;
            *.fas)   stem="\${stem%.fas}" ;;
        esac

        case "\$stem" in
            USERREF_*)
                echo "[REFERENCE_PREP] Reference genome stem '\$stem' (from '\$b') starts with the reserved prefix 'USERREF_'. Rename the reference file." >&2
                exit 1 ;;
        esac

        if grep -qxF "\$stem" seen_stems.txt; then
            echo "[REFERENCE_PREP] Duplicate reference genome stem '\$stem' (from '\$b'). Reference basenames must be unique." >&2
            exit 1
        fi
        echo "\$stem" >> seen_stems.txt

        out="reference_genomes/\${stem}.fasta"
        case "\$b" in
            *.gz) gzip -dc "\$f" > "\$out" ;;
            *)    cp -L "\$f" "\$out" ;;
        esac
        cp -L "\$out" "gtdbtk_refs/USERREF_\${stem}.fasta"
    done

    if ! ls reference_genomes/*.fasta >/dev/null 2>&1; then
        echo "[REFERENCE_PREP] No reference genome FASTAs were found in the supplied directory." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //;s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p reference_genomes gtdbtk_refs
    for f in refs_in/*; do
        [ -e "\$f" ] || continue
        b=\$(basename "\$f"); stem="\${b%.*}"; case "\$stem" in *.f*) stem="\${stem%.*}" ;; esac
        echo ">\${stem}_1" > "reference_genomes/\${stem}.fasta"
        echo "ACGT"       >> "reference_genomes/\${stem}.fasta"
    done
    ls reference_genomes/*.fasta >/dev/null 2>&1 || { echo ">ref_1" > reference_genomes/ref.fasta; echo "ACGT" >> reference_genomes/ref.fasta; }
    for g in reference_genomes/*.fasta; do
        cp -L "\$g" "gtdbtk_refs/USERREF_\$(basename "\$g")"
    done
    echo '"${task.process}": {bash: stub}' > versions.yml
    """
}

/*
 * Validate a user-supplied CheckM2 quality_report.tsv: every standardized
 * reference genome stem must appear in the report's Name column, otherwise
 * hard-fail listing the missing ones. Emits the validated report unchanged.
 */
process REFERENCE_CHECKM2_VALIDATE {
    label 'process_single'

    input:
    path(genomes, stageAs: 'reference_genomes/*')
    path(report)

    output:
    path 'reference_checkm2_report.tsv', emit: report
    path 'versions.yml',                 emit: versions

    script:
    """
    # Collect standardized reference stems.
    : > ref_stems.txt
    for f in reference_genomes/*.fasta; do
        [ -e "\$f" ] || continue
        b=\$(basename "\$f"); echo "\${b%.fasta}" >> ref_stems.txt
    done

    # Collect the Name column from the CheckM2 report.
    awk -F '\\t' 'NR==1 { for(i=1;i<=NF;i++) if(\$i=="Name"){c=i} next } (c){ print \$c }' "${report}" > report_names.txt
    if [ ! -s report_names.txt ]; then
        echo "[REFERENCE_CHECKM2_VALIDATE] Could not find a 'Name' column in ${report}." >&2
        exit 1
    fi

    # Any reference stem not present in the report is a hard error.
    missing=\$(grep -vxF -f report_names.txt ref_stems.txt || true)
    if [ -n "\$missing" ]; then
        echo "[REFERENCE_CHECKM2_VALIDATE] The supplied CheckM2 report (--reference_genomes_checkm2) is missing these reference genomes:" >&2
        echo "\$missing" | sed 's/^/  - /' >&2
        echo "Every reference genome must appear in the report's Name column. Add them, or omit --reference_genomes_checkm2 to run CheckM2 on the references." >&2
        exit 1
    fi

    cp -L "${report}" reference_checkm2_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1 | sed 's/.*[Aa]wk //;s/,.*//' || echo NA)
    END_VERSIONS
    """

    // No stub: this is a tool-free data-validation gate (bash/awk/grep only), so the
    // real check runs under -stub too — a bad --reference_genomes_checkm2 must fail early.
}
