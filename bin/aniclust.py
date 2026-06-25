#!/usr/bin/env python
"""
Greedy, UCLUST-like clustering of sequences from an ANI table.

Faithful reimplementation of CheckV's aniclust.py (BSD-licensed, Nayfach et al.),
vendored so the CHECKV_CLUSTER step is self-contained. Sequences are sorted by
length (longest first) and become cluster representatives; remaining sequences
join a representative if ANI/coverage thresholds are met (MIUVIG defaults:
--min_ani 95 --min_tcov 85 --min_qcov 0).
"""
import argparse
from collections import defaultdict


def parse_seqs(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        if name:
            yield name, ''.join(seq)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fna', required=True, help='FASTA of sequences that were compared')
    p.add_argument('--ani', required=True, help='ANI tsv from anicalc.py')
    p.add_argument('--out', required=True, help='output cluster tsv (representative<TAB>members)')
    p.add_argument('--min_ani', type=float, default=95.0)
    p.add_argument('--min_tcov', type=float, default=85.0)
    p.add_argument('--min_qcov', type=float, default=0.0)
    args = p.parse_args()

    lengths = {name: len(seq) for name, seq in parse_seqs(args.fna)}

    # edges[query] = set(targets passing thresholds)
    edges = defaultdict(set)
    with open(args.ani) as fh:
        next(fh)  # header
        for line in fh:
            q, t, _num, ani, qcov, tcov = line.rstrip('\n').split('\t')
            if q == t:
                continue
            if float(ani) >= args.min_ani and float(tcov) >= args.min_tcov and float(qcov) >= args.min_qcov:
                edges[q].add(t)

    clustered = set()
    with open(args.out, 'w') as out:
        for name in sorted(lengths, key=lambda n: lengths[n], reverse=True):
            if name in clustered:
                continue
            members = [name] + [t for t in edges[name] if t not in clustered]
            for m in members:
                clustered.add(m)
            out.write(name + '\t' + ','.join(members) + '\n')


if __name__ == '__main__':
    main()
