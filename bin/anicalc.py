#!/usr/bin/env python
"""
Compute pairwise ANI + alignment fraction from an all-vs-all blastn table.

Faithful reimplementation of CheckV's anicalc.py (BSD-licensed, Nayfach et al.),
vendored so the CHECKV_CLUSTER step is self-contained. Input is blastn outfmt
'6 std qlen slen'. Output columns: qname, tname, num_alns, ani, qcov, tcov.
"""
import argparse
from collections import defaultdict
from itertools import groupby


def parse_blast(handle):
    fields = ['qname', 'tname', 'pid', 'len', 'mm', 'gaps',
              'qstart', 'qstop', 'tstart', 'tstop', 'eval', 'score',
              'qlen', 'tlen']
    types = [str, str, float, int, int, int, int, int, int, int, float, float, int, int]
    for line in handle:
        values = line.split('\t')
        rec = dict((f, t(v)) for f, t, v in zip(fields, types, values))
        yield rec


def yield_alignment_blocks(handle):
    # group consecutive lines by (query, target) pair
    key = lambda r: (r['qname'], r['tname'])
    for _, alns in groupby(parse_blast(handle), key=key):
        yield list(alns)


def compute_ani(alns):
    return round(sum(a['len'] * a['pid'] for a in alns) / sum(a['len'] for a in alns), 2)


def compute_cov(alns, which):
    # compute the fraction of the (q or t) sequence covered by alignments,
    # merging overlapping intervals.
    coords = sorted([(a['%sstart' % which], a['%sstop' % which]) for a in alns])
    nr = [list(coords[0])]
    for start, stop in coords[1:]:
        if start <= nr[-1][1]:
            nr[-1][1] = max(nr[-1][1], stop)
        else:
            nr.append([start, stop])
    alen = sum(stop - start + 1 for start, stop in nr)
    seqlen = alns[0]['%slen' % which]
    return round(100.0 * alen / seqlen, 2)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i', dest='blast', required=True, help="blastn tsv (outfmt '6 std qlen slen')")
    p.add_argument('-o', dest='out', required=True, help='output ANI tsv')
    args = p.parse_args()

    with open(args.blast) as fh, open(args.out, 'w') as out:
        out.write('\t'.join(['qname', 'tname', 'num_alns', 'ani', 'qcov', 'tcov']) + '\n')
        for alns in yield_alignment_blocks(fh):
            qname, tname = alns[0]['qname'], alns[0]['tname']
            ani = compute_ani(alns)
            qcov = compute_cov(alns, 'q')
            tcov = compute_cov(alns, 't')
            out.write('\t'.join(map(str, [qname, tname, len(alns), ani, qcov, tcov])) + '\n')


if __name__ == '__main__':
    main()
