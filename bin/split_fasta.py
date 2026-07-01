#!/usr/bin/env python3
"""
Split a FASTA file into N chunks of roughly equal sequence count, preserving each
record byte-for-byte. Used to scatter a large gene catalogue across parallel DRAM
annotation tasks (the ORF annotations are independent, so the pieces can be merged
back afterwards).

Writes <prefix>.001.faa, <prefix>.002.faa, … (1-based, zero-padded).

Usage: split_fasta.py INPUT.faa N_PARTS PREFIX
"""
import math
import sys


def iter_records(path):
    """Yield (header_line, body_str) with original newlines preserved."""
    header = None
    body = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(body)
                header = line
                body = []
            else:
                body.append(line)
    if header is not None:
        yield header, ''.join(body)


def main():
    if len(sys.argv) != 4:
        sys.stderr.write(__doc__)
        return 2
    inp, n_parts, prefix = sys.argv[1], int(sys.argv[2]), sys.argv[3]

    records = list(iter_records(inp))
    total = len(records)
    if total == 0:
        sys.stderr.write("split_fasta: no sequences in %s\n" % inp)
        return 1

    n_parts = max(1, min(n_parts, total))
    size = math.ceil(total / n_parts)

    written = 0
    for start in range(0, total, size):
        written += 1
        chunk = records[start:start + size]
        with open("%s.%03d.faa" % (prefix, written), 'w') as out:
            for header, body in chunk:
                out.write(header)
                out.write(body)

    sys.stderr.write("split_fasta: %d sequences -> %d chunks (~%d each)\n"
                     % (total, written, size))
    return 0


if __name__ == '__main__':
    sys.exit(main())
