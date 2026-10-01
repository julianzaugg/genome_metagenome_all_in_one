#!/usr/bin/env bash
# Print the IDs of high-quality genomes: completeness - WEIGHT*contamination >= QUALITY.
#
# usage: select_hq_ids.sh SOURCE QUALITY WEIGHT CHECKM2_REPORT CHECKM1_REPORT
#   SOURCE   checkm2 | checkm1 | either (pass in at least one report)
#            | both (pass in every report; both must be non-empty)
#   reports  CheckM2 quality_report.tsv / CheckM1 tab table; an empty or missing
#            path means that report was not produced.
# Column positions are found by header name (Name / Bin Id / genome,
# Completeness, Contamination), so either tool's layout works.
set -euo pipefail

source=$1; quality=$2; weight=$3; checkm2=${4:-}; checkm1=${5:-}

pass_ids() {
    [ -n "$1" ] && [ -s "$1" ] || return 0
    awk -F '\t' -v q="$quality" -v k="$weight" '
        NR==1 { for (i=1; i<=NF; i++) { if ($i=="Completeness") cc=i; if ($i=="Contamination") ct=i;
                                        if ($i=="Name" || $i=="Bin Id" || $i=="genome") id=i } next }
        (cc && ct && id && ($cc - k*$ct) >= q) { print $id }
    ' "$1" | sort -u
}

case "$source" in
    checkm2) pass_ids "$checkm2" ;;
    checkm1) pass_ids "$checkm1" ;;
    either)  { pass_ids "$checkm2"; pass_ids "$checkm1"; } | sort -u ;;
    both)
        if [ ! -s "$checkm2" ] || [ ! -s "$checkm1" ]; then
            echo "select_hq_ids.sh: source 'both' needs both a CheckM2 and a CheckM1 report" >&2
            exit 1
        fi
        comm -12 <(pass_ids "$checkm2") <(pass_ids "$checkm1")
        ;;
    *) echo "select_hq_ids.sh: unknown source '$source' (checkm1|checkm2|either|both)" >&2; exit 1 ;;
esac
