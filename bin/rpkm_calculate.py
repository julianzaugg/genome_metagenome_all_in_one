#!/usr/bin/env python3
"""Calculate gene-catalogue RPKM normalized by SingleM marker-gene RPKM."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row:
                yield row


def write_tsv(path, fieldnames, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def as_float(value, default=0.0):
    try:
        if value is None or value == "":
            return default
        number = float(str(value).replace(",", ""))
        if math.isnan(number) or math.isinf(number):
            return default
        return number
    except ValueError:
        return default


def load_marker_lengths(path):
    lengths = {}
    for row in read_tsv(path):
        gene = row.get("Gene")
        if gene:
            lengths[gene] = round(as_float(row.get("avg_len")))
    return lengths


def marker_name_from_blast(path):
    name = Path(path).name
    return name[:-10] if name.endswith("_blast.tsv") else Path(path).stem


def calculate_singlem(singlem_blast_dir, marker_lengths):
    sample_marker_counts = defaultdict(int)
    sample_totals = defaultdict(int)

    for blast_path in sorted(Path(singlem_blast_dir).rglob("*_blast.tsv")):
        marker = marker_name_from_blast(blast_path)
        for row in read_tsv(blast_path):
            sample = row.get("sample")
            if not sample:
                continue
            sample_marker_counts[(sample, marker)] += 1
            sample_totals[sample] += 1

    singlem_rows = []
    sample_values = defaultdict(list)
    for sample, marker in sorted(sample_marker_counts):
        read_count = sample_marker_counts[(sample, marker)]
        total = sample_totals[sample]
        avg_len = marker_lengths.get(marker, 0)
        denominator = (avg_len * 3 / 1000) * (total / 1_000_000)
        rpkm = read_count / denominator if denominator else 0.0
        singlem_rows.append({"sample": sample, "Marker_gene": marker, "rpkm": f"{rpkm:.12g}"})
        sample_values[sample].append(rpkm)

    mean_rows = []
    means = {}
    for sample in sorted(sample_values):
        values = sample_values[sample]
        mean = sum(values) / len(values) if values else 0.0
        means[sample] = mean
        mean_rows.append({"sample": sample, "Mean_rpkm": f"{mean:.12g}"})

    return singlem_rows, mean_rows, means


def calculate_gene_tables(gene_blast, singlem_means):
    gene_counts = defaultdict(int)
    gene_lengths = {}
    sample_totals = defaultdict(int)
    samples = set()
    genes = set()

    for row in read_tsv(gene_blast):
        sample = row.get("sample")
        gene = row.get("sseqid")
        if not sample or not gene:
            continue
        length = as_float(row.get("slen"))
        key = (gene, sample)
        gene_counts[key] += 1
        gene_lengths.setdefault(gene, length)
        sample_totals[sample] += 1
        samples.add(sample)
        genes.add(gene)

    ordered_samples = sorted(samples)
    count_rows = []
    norm_rows = []

    for gene in sorted(genes):
        count_row = {"Gene_ID": gene}
        norm_row = {"Gene_ID": gene}
        gene_len = gene_lengths.get(gene, 0.0)

        for sample in ordered_samples:
            count = gene_counts.get((gene, sample), 0)
            count_row[sample] = str(count)

            total = sample_totals.get(sample, 0)
            denominator = (gene_len * 3 / 1000) * (total / 1_000_000)
            rpkm = count / denominator if denominator else 0.0
            singlem_mean = singlem_means.get(sample, 0.0)
            normalised = (rpkm / singlem_mean) * 100 if singlem_mean else 0.0
            norm_row[sample] = f"{normalised:.12g}"

        count_rows.append(count_row)
        norm_rows.append(norm_row)

    fields = ["Gene_ID"] + ordered_samples
    return fields, norm_rows, count_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--singlem-blast-dir", required=True)
    parser.add_argument("--marker-lengths", required=True)
    parser.add_argument("--gene-blast", required=True)
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    marker_lengths = load_marker_lengths(args.marker_lengths)
    singlem_rows, mean_rows, singlem_means = calculate_singlem(args.singlem_blast_dir, marker_lengths)
    gene_fields, norm_rows, count_rows = calculate_gene_tables(args.gene_blast, singlem_means)

    write_tsv(outdir / "singlem_sample_rpkm.tsv", ["sample", "Marker_gene", "rpkm"], singlem_rows)
    write_tsv(outdir / "singlem_rpkm_means.tsv", ["sample", "Mean_rpkm"], mean_rows)
    write_tsv(outdir / "gene_catalogue_rpkm_per_gene_normalised.tsv", gene_fields, norm_rows)
    write_tsv(outdir / "gene_catalogue_mapped_reads_per_gene.tsv", gene_fields, count_rows)


if __name__ == "__main__":
    main()
