#!/usr/bin/env python3
"""
Reshape `inStrain compare` output into per-sample-pair strain-sharing calls.

`inStrain compare` answers, for every genome shared by two samples, how similar the
populations mapping to that genome are. The two numbers that matter are:

  popANI                   population-level ANI. Only counts a position as a difference
                           when the two samples share NO alleles there -- so a position
                           where sample A is fixed A and sample B is polymorphic A/G is
                           NOT a difference. This is what makes popANI a strain-sharing
                           statistic rather than a consensus-similarity one.
  percent_genome_compared  fraction of the genome covered well enough in BOTH samples to
                           be compared at all. popANI over a tiny fraction of a genome is
                           not evidence of anything, which is why inStrain gates on it.

The conventional same-strain call is popANI >= 0.99999 AND percent_genome_compared >= 0.5.
This script does NOT re-apply those thresholds: `inStrain compare` already clustered the
genomes under whatever -ani/-cov it was given (see conf/modules.config), and its
strain_clusters output is the authority. Re-deriving the call here would create a second
source of truth that silently disagrees whenever the ext.args are tuned. Two samples are
reported as sharing a strain iff they land in the same cluster for that genome.

Outputs (under --outdir):

  strain_sharing_summary.tsv   one row per (genome, sample_a, sample_b): popANI, conANI,
                               percent_genome_compared, coverage_overlap, cluster ids and
                               the same_strain call
  strain_sharing_counts.tsv    one row per sample pair: how many genomes were comparable
                               and how many of those shared a strain. This is the table to
                               read first -- it is the "do these two samples share
                               populations" answer, aggregated over all genomes.
  strain_sharing_matrix/<genome>.popani.tsv
                               sample x sample popANI matrix per genome, for plotting

A caveat worth carrying into interpretation: a pair only appears here for genomes with
enough coverage in BOTH samples. Absence of a row is "not enough data", not "different
strains" -- the two are easy to conflate when scanning these tables.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def find_one(directory, suffix):
    """Locate a single inStrain output file by suffix, searching output/ then the root."""
    for sub in ("output", ""):
        d = os.path.join(directory, sub) if sub else directory
        if not os.path.isdir(d):
            continue
        hits = sorted(f for f in os.listdir(d) if f.endswith(suffix))
        if hits:
            return os.path.join(d, hits[0])
    return ""


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def get(row, *names, default=""):
    """inStrain column names have drifted across versions; accept any of `names`."""
    for n in names:
        if n in row and row[n] != "":
            return row[n]
    return default


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def pair_key(a, b):
    return tuple(sorted((a, b)))


def load_clusters(path):
    """genome -> {sample: cluster_id} from inStrain's strain_clusters table."""
    clusters = defaultdict(dict)
    for row in read_tsv(path):
        genome = get(row, "genome")
        sample = get(row, "sample", "name")
        cluster = get(row, "cluster", "strain_cluster")
        if genome and sample:
            clusters[genome][sample] = cluster
    return clusters


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--compare-dir", required=True,
                    help="inStrain compare output directory (the *.IS dir)")
    ap.add_argument("--outdir", default="summary", help="output directory")
    ap.add_argument("--genome-wide", default="",
                    help="explicit path to *_genomeWide_compare.tsv (default: found under --compare-dir)")
    ap.add_argument("--strain-clusters", default="",
                    help="explicit path to *_strain_clusters.tsv (default: found under --compare-dir)")
    args = ap.parse_args()

    gw_path = args.genome_wide or find_one(args.compare_dir, "genomeWide_compare.tsv")
    sc_path = args.strain_clusters or find_one(args.compare_dir, "strain_clusters.tsv")

    if not gw_path:
        sys.exit(f"[instrain_summarise] no *_genomeWide_compare.tsv found under {args.compare_dir}")
    if not sc_path:
        # Clustering needs both an .stb and -ani; without it we can still report the
        # numbers, just not inStrain's own same-strain verdict.
        print(f"[instrain_summarise] WARNING: no *_strain_clusters.tsv under {args.compare_dir}; "
              "same_strain will be NA", file=sys.stderr)

    clusters = load_clusters(sc_path)
    rows = read_tsv(gw_path)

    os.makedirs(args.outdir, exist_ok=True)
    matrix_dir = os.path.join(args.outdir, "strain_sharing_matrix")
    os.makedirs(matrix_dir, exist_ok=True)

    summary_path = os.path.join(args.outdir, "strain_sharing_summary.tsv")
    counts = defaultdict(lambda: {"compared": 0, "same": 0})
    popani = defaultdict(dict)   # genome -> {(a,b): popANI}
    samples = defaultdict(set)   # genome -> {sample}

    with open(summary_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "genome", "sample_a", "sample_b", "popANI", "conANI",
            "percent_genome_compared", "coverage_overlap",
            "cluster_a", "cluster_b", "same_strain",
        ])

        for row in rows:
            genome = get(row, "genome")
            a = get(row, "name1", "sample1")
            b = get(row, "name2", "sample2")
            if not genome or not a or not b or a == b:
                continue

            ca = clusters.get(genome, {}).get(a, "")
            cb = clusters.get(genome, {}).get(b, "")
            if ca and cb:
                same = "TRUE" if ca == cb else "FALSE"
            else:
                same = "NA"

            writer.writerow([
                genome, a, b,
                get(row, "popANI"), get(row, "conANI"),
                get(row, "percent_genome_compared"),
                get(row, "coverage_overlap"),
                ca or "NA", cb or "NA", same,
            ])

            key = pair_key(a, b)
            counts[key]["compared"] += 1
            if same == "TRUE":
                counts[key]["same"] += 1

            value = to_float(get(row, "popANI"))
            if value is not None:
                popani[genome][key] = value
            samples[genome].update((a, b))

    counts_path = os.path.join(args.outdir, "strain_sharing_counts.tsv")
    with open(counts_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_a", "sample_b", "n_genomes_compared", "n_same_strain"])
        for (a, b), c in sorted(counts.items()):
            writer.writerow([a, b, c["compared"], c["same"]])

    for genome, pairs in popani.items():
        names = sorted(samples[genome])
        # '/' and os.sep would fan the matrices out into unintended subdirectories.
        safe = genome.replace(os.sep, "_").replace("/", "_")
        path = os.path.join(matrix_dir, f"{safe}.popani.tsv")
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample"] + names)
            for a in names:
                line = [a]
                for b in names:
                    if a == b:
                        line.append("1.0")
                    else:
                        v = pairs.get(pair_key(a, b))
                        line.append("NA" if v is None else f"{v:.7f}")
                writer.writerow(line)

    print(f"[instrain_summarise] {len(rows)} comparisons, {len(counts)} sample pairs, "
          f"{len(popani)} genomes -> {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
