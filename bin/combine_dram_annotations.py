#!/usr/bin/env python3
"""
Concatenate multiple DRAM annotations.tsv files into one (row-wise, union of columns).

DRAM writes annotations.tsv with the gene id as the (unnamed) first/index column and
a per-genome 'fasta' column. Different inputs may carry different columns — a database
with no hits for a given genome contributes no column — so we outer-align on the union
of columns and leave missing cells empty.

Two uses in this pipeline:
  * merge per-bin annotations into one table so a single DRAM distill produces a
    cross-MAG summary (distill groups rows by the 'fasta' column); and
  * reassemble the pieces of a gene catalogue that was split for parallel annotation.

Gene ids are unique across inputs (DRAM prefixes them per fasta / they are catalogue
rep ids), so the row index stays unique on concatenation.

Usage:
  combine_dram_annotations.py OUT.tsv IN1.tsv IN2.tsv ...
  combine_dram_annotations.py OUT.tsv --glob 'anno_*.tsv'
"""
import glob
import sys

import pandas as pd


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        return 2
    out = sys.argv[1]
    rest = sys.argv[2:]
    if rest[0] == '--glob':
        inputs = sorted(glob.glob(rest[1]))
    else:
        inputs = rest
    if not inputs:
        sys.stderr.write("combine_dram_annotations: no input files\n")
        return 2

    # dtype=str keeps every value exactly as DRAM wrote it (no float coercion of
    # score columns) and makes the outer concat fill gaps with empty strings.
    frames = [pd.read_csv(f, sep='\t', index_col=0, dtype=str) for f in inputs]
    combined = pd.concat(frames, axis=0, sort=False)

    # Keep 'fasta' as the leading column, as DRAM emits it, so distill's groupby
    # and any downstream reader see the expected layout.
    if 'fasta' in combined.columns:
        cols = ['fasta'] + [c for c in combined.columns if c != 'fasta']
        combined = combined[cols]

    combined.to_csv(out, sep='\t')
    sys.stderr.write(
        "combine_dram_annotations: %d files -> %s (%d genes, %d columns)\n"
        % (len(inputs), out, combined.shape[0], combined.shape[1]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
