#!/usr/bin/env python3
"""
Assess how well a reference (assembly / MAG set) represents sequencing data,
from BAM files alone.

For short reads, "% of reads mapped" is a fine summary. For long reads it is
misleading, because reads vary enormously in length, so a read count says little
about how much *sequence data* was represented. The base-level metric
`mapped read bases / total bases sequenced` is what you want there.

This tool reports three base-level metrics, per (sample, reference set):

  reads_mapped          number of primary mapped reads
  pct_reads_mapped      reads_mapped / total reads sequenced

  bases_mapped_readlen  metric A: full SEQ length of every primary mapped read
                        (each mapped read counted once at its full length --
                        "how many bases of read data mapped somewhere")
  pct_bases_mapped_A    bases_mapped_readlen / total bases sequenced

  bases_mapped_cigar    metric B: aligned bases only (CIGAR M/I/=/X) of PRIMARY
                        alignments, excluding soft-clipped/unmapped tails AND
                        excluding supplementary (split) alignment segments
  pct_bases_mapped_B    bases_mapped_cigar / total bases sequenced

  bases_mapped_cigar_all metric C: aligned bases (CIGAR M/I/=/X) of primary AND
                        supplementary alignments -- the honest aligned-base
                        total for long reads that split across a contig/genome
                        boundary or a chimera. reads_supplementary counts the
                        supplementary records contributing to it.
  pct_bases_mapped_C    bases_mapped_cigar_all / total bases sequenced

Provenance: ported from https://github.com/julianzaugg/assess_long_read_mapping
(metrics A and B), with metric C added here. That repo's own doc claims
`samtools stats`' "bases mapped (cigar)" line is primary-alignment-only; this is
NOT reliable in general -- empirically (samtools 1.22-1.24, real CoverM/minimap2
BAMs) it silently includes supplementary CIGAR bases whenever the supplementary
record carries a real (hard-clipped) SEQ, which is minimap2's default. So B here
is always computed with an explicit `-F 0x900` (or from a pre-filtered
--stats-file/--per-contig input), never from a bare `samtools stats` call.

Two ways to get the whole-BAM SN numbers (reads_mapped, reads_total,
reads_unmapped, bases_total, bases_mapped, bases_mapped (cigar)):
  * pass a BAM (positional / --bam-dir): this script shells out to
    `samtools stats -F 0x900 <bam>` itself.
  * pass --stats-file: consume already-computed `samtools stats -F 0x900`
    output (the pipeline runs samtools once and passes the text through,
    rather than paying for it twice).

Metric C and the per-genome breakdown come from a per-contig CIGAR table
(--per-contig): one row per contig with reads_primary, reads_supplementary,
bases_A, bases_B, bases_C, produced by a `samtools view -F 0x104 | awk` pass
over primary+supplementary records (excluding secondary and unmapped),
parsing CIGAR directly (never SEQ length) so it is correct regardless of
aligner SEQ conventions. When running against a bare BAM instead (no
--per-contig given), this script computes the same figures itself in pure
Python (stdlib only) via `samtools view -F 0x104 <bam>`.

The denominators (total reads / total bases sequenced) are handled automatically:
  * If the BAM still contains unmapped reads, they are taken from the BAM itself
    (`raw total sequences` / `total length`) and total_source = "bam".
  * Otherwise (e.g. a `coverm filter` BAM, or CoverM's own `--discard-unmapped`,
    both of which drop unmapped reads) they are taken from a supplied --totals
    table, or --totals-from-seqkit (a `seqkit stats --tabular --all` file --
    summed across rows, e.g. R1+R2), and total_source = "supplied". With
    neither, the tool falls back to the BAM's own totals (yielding ~100%) and
    warns that the percentage is relative only to the reads present in the BAM.

The tool never reproduces CoverM's coverage maths; it only reads properties of
the BAM. Run it against a `coverm filter` BAM and the alignments are already
identity/aligned-% filtered, so the numbers reflect that filtering. Run it
against a raw/cached BAM and every primary alignment is counted (no identity
filter), so % mapped will read higher -- choose the input BAM deliberately.
"""
import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys

# samtools stats SN keys -> the field name we use internally.
SN_KEYS = {
    "raw total sequences": "reads_total",
    "reads mapped": "reads_mapped",
    "reads unmapped": "reads_unmapped",
    "bases mapped": "bases_mapped_readlen",       # metric A (full read length)
    "bases mapped (cigar)": "bases_mapped_cigar",  # metric B, IF -F 0x900 was used
    "total length": "bases_total",
}

# Flags samtools stats must be run with so "bases mapped (cigar)" is strictly
# primary-only (see the provenance note above for why this cannot be assumed
# from a bare `samtools stats` call).
STRICT_PRIMARY_FLAG = "0x900"

OUT_COLUMNS = [
    "sample",
    "set",
    "reads_mapped",
    "reads_total",
    "pct_reads_mapped",
    "reads_supplementary",
    "bases_mapped_readlen",
    "pct_bases_mapped_A",
    "bases_mapped_cigar",
    "pct_bases_mapped_B",
    "bases_mapped_cigar_all",
    "pct_bases_mapped_C",
    "bases_total",
    "total_source",
    "note",
]

PER_GENOME_COLUMNS = [
    "set",
    "sample",
    "genome",
    "reads_primary",
    "reads_supplementary",
    "bases_A",
    "pct_bases_mapped_A",
    "bases_B",
    "pct_bases_mapped_B",
    "bases_C",
    "pct_bases_mapped_C",
]

CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")
ALIGNED_OPS = set("MI=X")


def warn(msg):
    print(f"[assess_mapping] WARNING: {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# samtools stats (whole-BAM A / B / denominators)
# ---------------------------------------------------------------------------
def run_samtools_stats(samtools, bam, threads):
    """Run `samtools stats -F 0x900` (strict primary-only) and return SN fields."""
    cmd = [samtools, "stats", "-F", STRICT_PRIMARY_FLAG]
    if threads and threads > 1:
        cmd += ["-@", str(threads)]
    cmd.append(bam)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"samtools stats failed on {bam}")
    return parse_sn(proc.stdout)


def parse_sn(text):
    """
    Parse the SN (summary numbers) block of samtools stats output.

    Lines look like:  SN\t<key>:\t<value>\t# optional comment
    """
    fields = {}
    for line in text.splitlines():
        if not line.startswith("SN\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        key = parts[1].rstrip(":").strip()
        if key not in SN_KEYS:
            continue
        try:
            fields[SN_KEYS[key]] = int(float(parts[2]))
        except ValueError:
            continue
    return fields


# ---------------------------------------------------------------------------
# Metric C / per-contig CIGAR pass (pure Python fallback for standalone BAM use)
# ---------------------------------------------------------------------------
def cigar_aligned_bases(cigar):
    """Sum of CIGAR M/I/=/X op lengths -- the read bases that took part in an
    alignment. Parses the CIGAR string directly, never SEQ length, so it is
    correct whether or not the aligner wrote a real SEQ for this record."""
    if cigar in ("", "*"):
        return 0
    total = 0
    for length, op in CIGAR_OP_RE.findall(cigar):
        if op in ALIGNED_OPS:
            total += int(length)
    return total


def compute_per_contig_from_bam(samtools, bam, threads):
    """
    Fallback for standalone BAM use (no --per-contig supplied): reproduce the
    pipeline's own `samtools view -F 0x104` (exclude secondary + unmapped, keep
    supplementary) CIGAR pass in pure Python. Returns {contig: {'reads_primary',
    'reads_supplementary', 'bases_A', 'bases_B', 'bases_C'}}.
    """
    cmd = [samtools, "view", "-F", "0x104"]
    if threads and threads > 1:
        cmd += ["-@", str(threads)]
    cmd.append(bam)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or f"samtools view failed on {bam}")

    per_contig = {}
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        rname = fields[2]
        cigar = fields[5]
        seq = fields[9]
        is_supplementary = bool(flag & 0x800)
        bases = cigar_aligned_bases(cigar)
        rec = per_contig.setdefault(rname, {
            "reads_primary": 0, "reads_supplementary": 0,
            "bases_A": 0, "bases_B": 0, "bases_C": 0,
        })
        if is_supplementary:
            rec["reads_supplementary"] += 1
            rec["bases_C"] += bases
        else:
            rec["reads_primary"] += 1
            rec["bases_A"] += 0 if seq == "*" else len(seq)
            rec["bases_B"] += bases
            rec["bases_C"] += bases
    return per_contig


def read_per_contig_file(path):
    """Read a precomputed per-contig CIGAR TSV (contig, reads_primary,
    reads_supplementary, bases_A, bases_B, bases_C) -> {contig: {...}}."""
    per_contig = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            per_contig[row["contig"]] = {
                "reads_primary": int(row.get("reads_primary") or 0),
                "reads_supplementary": int(row.get("reads_supplementary") or 0),
                "bases_A": int(row.get("bases_A") or 0),
                "bases_B": int(row.get("bases_B") or 0),
                "bases_C": int(row.get("bases_C") or 0),
            }
    return per_contig


def read_contig_map(path):
    """Read a two-column (no header) contig -> genome TSV."""
    contig_genome = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            contig_genome[parts[0]] = parts[1]
    return contig_genome


def aggregate_per_genome(per_contig, contig_genome):
    """Group per-contig rows by genome (via contig_genome). Contigs absent
    from the map are pooled under 'unassigned' and a warning is issued."""
    per_genome = {}
    unmapped_contigs = []
    for contig, rec in per_contig.items():
        genome = contig_genome.get(contig)
        if genome is None:
            genome = "unassigned"
            unmapped_contigs.append(contig)
        g = per_genome.setdefault(genome, {
            "reads_primary": 0, "reads_supplementary": 0,
            "bases_A": 0, "bases_B": 0, "bases_C": 0,
        })
        for k in ("reads_primary", "reads_supplementary", "bases_A", "bases_B", "bases_C"):
            g[k] += rec[k]
    if unmapped_contigs:
        warn(f"{len(unmapped_contigs)} contig(s) not found in --contig-map, "
             f"pooled as 'unassigned': {', '.join(sorted(unmapped_contigs)[:5])}"
             + (" ..." if len(unmapped_contigs) > 5 else ""))
    return per_genome


# ---------------------------------------------------------------------------
# Denominator sources
# ---------------------------------------------------------------------------
def sample_name(bam, strip_suffixes):
    base = os.path.basename(bam)
    if base.endswith(".bam"):
        base = base[: -len(".bam")]
    for suf in strip_suffixes:
        if suf and base.endswith(suf):
            base = base[: -len(suf)]
    return base


def read_totals(path):
    """Read a totals TSV -> {sample: {'total_reads': int, 'total_bases': int}}."""
    totals = {}
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    if not rows:
        return totals
    header = [h.strip().lower() for h in rows[0]]

    def idx(*names):
        for n in names:
            if n in header:
                return header.index(n)
        return None

    s_i = idx("sample", "sample_id", "id")
    r_i = idx("total_reads", "reads", "num_seqs", "raw_count")
    b_i = idx("total_bases", "bases", "sum_len", "bp")
    if s_i is None:
        sys.exit(f"--totals: no 'sample' column found in {path}")
    for row in rows[1:]:
        if s_i >= len(row):
            continue
        rec = {}
        if r_i is not None and r_i < len(row) and row[r_i].strip():
            rec["total_reads"] = int(float(row[r_i]))
        if b_i is not None and b_i < len(row) and row[b_i].strip():
            rec["total_bases"] = int(float(row[b_i]))
        totals[row[s_i].strip()] = rec
    return totals


def read_totals_from_seqkit(path):
    """
    Read a `seqkit stats --tabular --all` TSV and sum num_seqs/sum_len across
    ALL rows (paired reads => one row per mate file) -> a single
    {'total_reads': int, 'total_bases': int} record for this sample.
    """
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    if not rows:
        return {}
    header = rows[0]

    def col_index(suffix):
        for i, name in enumerate(header):
            if name == suffix or name.endswith(" " + suffix):
                return i
        return None

    n_i = col_index("num_seqs")
    b_i = col_index("sum_len")
    reads = bases = 0
    for row in rows[1:]:
        if n_i is not None and n_i < len(row):
            reads += int(float(row[n_i]))
        if b_i is not None and b_i < len(row):
            bases += int(float(row[b_i]))
    return {"total_reads": reads, "total_bases": bases}


def collect_bams(inputs, bam_dir):
    """Expand positional inputs (files/globs/dirs) plus --bam-dir into BAM paths."""
    bams = []
    candidates = list(inputs)
    if bam_dir:
        candidates.append(bam_dir)
    for item in candidates:
        if os.path.isdir(item):
            bams += sorted(glob.glob(os.path.join(item, "*.bam")))
        elif any(ch in item for ch in "*?["):
            bams += sorted(glob.glob(item))
        else:
            bams.append(item)
    # de-duplicate, preserve order
    seen = set()
    out = []
    for b in bams:
        real = os.path.realpath(b)
        if real not in seen:
            seen.add(real)
            out.append(b)
    return out


def pct(num, den):
    if not den:
        return ""
    return f"{num / den * 100:.2f}"


def resolve_denominators(fields, supplied):
    """
    Resolve (reads_total, bases_total, source, note) from the whole-BAM SN
    fields and an optional supplied totals record.
    """
    unmapped = fields.get("reads_unmapped", 0)
    bam_reads_total = fields.get("reads_total", 0)
    bam_bases_total = fields.get("bases_total", 0)

    if unmapped > 0:
        # BAM retains unmapped reads -> it is a complete record of the sample.
        return bam_reads_total, bam_bases_total, "bam", ""
    if supplied and ("total_reads" in supplied or "total_bases" in supplied):
        reads_total = supplied.get("total_reads", bam_reads_total)
        bases_total = supplied.get("total_bases", bam_bases_total)
        return reads_total, bases_total, "supplied", ""
    note = ("no unmapped reads in BAM and no --totals/--totals-from-seqkit; % is "
            "relative to reads in the BAM, not all sequenced data")
    return bam_reads_total, bam_bases_total, "bam(mapped-only)", note


def assess_one(fields, supplied, per_contig=None):
    """
    Given parsed SN fields, an optional supplied totals record, and an optional
    per-contig table (for metric C / reads_supplementary), resolve the
    denominators and compute the output row values. Returns (values, note).
    """
    reads_mapped = fields.get("reads_mapped", 0)
    a_bases = fields.get("bases_mapped_readlen", 0)
    b_bases = fields.get("bases_mapped_cigar", 0)

    if per_contig:
        c_bases = sum(r["bases_C"] for r in per_contig.values())
        reads_supp = sum(r["reads_supplementary"] for r in per_contig.values())
        sum_a = sum(r["bases_A"] for r in per_contig.values())
        sum_b = sum(r["bases_B"] for r in per_contig.values())
        if sum_a != a_bases:
            warn(f"self-check failed: per-contig bases_A sum ({sum_a}) != "
                 f"samtools stats bases_mapped ({a_bases})")
        if sum_b != b_bases:
            warn(f"self-check failed: per-contig bases_B sum ({sum_b}) != "
                 f"samtools stats bases_mapped_cigar ({b_bases})")
    else:
        # No per-genome breakdown available: C degrades to B (no visibility
        # into supplementary alignments without a CIGAR pass).
        c_bases = b_bases
        reads_supp = 0

    reads_total, bases_total, source, note = resolve_denominators(fields, supplied)
    if source == "bam(mapped-only)":
        warn(note)

    values = {
        "reads_mapped": reads_mapped,
        "reads_total": reads_total,
        "pct_reads_mapped": pct(reads_mapped, reads_total),
        "reads_supplementary": reads_supp,
        "bases_mapped_readlen": a_bases,
        "pct_bases_mapped_A": pct(a_bases, bases_total),
        "bases_mapped_cigar": b_bases,
        "pct_bases_mapped_B": pct(b_bases, bases_total),
        "bases_mapped_cigar_all": c_bases,
        "pct_bases_mapped_C": pct(c_bases, bases_total),
        "bases_total": bases_total,
        "total_source": source,
        "note": note,
    }
    return values, bases_total


def write_per_genome(path, per_genome, sample, set_label, bases_total):
    lines = ["\t".join(PER_GENOME_COLUMNS)]
    for genome in sorted(per_genome):
        rec = per_genome[genome]
        row = [
            set_label, sample, genome,
            str(rec["reads_primary"]), str(rec["reads_supplementary"]),
            str(rec["bases_A"]), pct(rec["bases_A"], bases_total),
            str(rec["bases_B"]), pct(rec["bases_B"], bases_total),
            str(rec["bases_C"]), pct(rec["bases_C"], bases_total),
        ]
        lines.append("\t".join(row))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("bams", nargs="*", help="BAM files, globs, or directories")
    ap.add_argument("--bam-dir", default="", help="directory of *.bam to include")
    ap.add_argument("--totals", default="",
                    help="TSV with columns sample, total_reads, total_bases "
                         "(the denominator source for filtered/unmapped-free BAMs)")
    ap.add_argument("--totals-from-seqkit", default="",
                    help="`seqkit stats --tabular --all` TSV; num_seqs/sum_len "
                         "summed across all rows becomes the denominator for "
                         "this one sample (pipeline mode)")
    ap.add_argument("--strip-suffix", action="append", default=[],
                    help="extra suffix to strip from the BAM basename to form the "
                         "sample id (repeatable), e.g. --strip-suffix .fastq.gz")
    ap.add_argument("-@", "--threads", type=int, default=1,
                    help="threads passed to samtools (default 1)")
    ap.add_argument("--samtools", default="samtools",
                    help="path to the samtools executable")
    ap.add_argument("-o", "--out", default="",
                    help="output TSV path (default: stdout)")

    # Pipeline-mode options: process exactly one (sample, set) combination.
    ap.add_argument("--sample", default="",
                    help="explicit sample id (pipeline mode; overrides the "
                         "BAM-filename-derived id)")
    ap.add_argument("--set", dest="set_label", default="",
                    help="reference-set label, added as the 'set' column "
                         "(e.g. HQ_Derep_MAGs)")
    ap.add_argument("--stats-file", default="",
                    help="pre-computed `samtools stats -F 0x900 <bam>` output; "
                         "when given, no BAM is needed for A/B/denominators")
    ap.add_argument("--per-contig", default="",
                    help="precomputed per-contig CIGAR TSV (contig, "
                         "reads_primary, reads_supplementary, bases_A, "
                         "bases_B, bases_C) for metric C + per-genome output")
    ap.add_argument("--contig-map", default="",
                    help="two-column (no header) contig -> genome TSV")
    ap.add_argument("--per-genome-out", default="",
                    help="output path for the per-genome breakdown table "
                         "(requires --per-contig and --contig-map)")
    args = ap.parse_args()

    samtools = shutil.which(args.samtools) or args.samtools
    pipeline_mode = bool(args.stats_file or args.sample or args.set_label)

    if not pipeline_mode and not shutil.which(samtools):
        sys.exit(f"samtools not found on PATH (looked for '{args.samtools}'). "
                 f"Activate the conda env that provides it, e.g. coverm_0.8.0.")

    supplied_totals = read_totals(args.totals) if args.totals else {}
    supplied_seqkit = (read_totals_from_seqkit(args.totals_from_seqkit)
                        if args.totals_from_seqkit else None)

    if pipeline_mode:
        if not args.stats_file:
            sys.exit("pipeline mode (--sample/--set) requires --stats-file")
        if not args.sample or not args.set_label:
            sys.exit("--stats-file requires both --sample and --set")
        with open(args.stats_file) as fh:
            fields = parse_sn(fh.read())
        if "reads_total" not in fields:
            sys.exit(f"no summary numbers found in --stats-file {args.stats_file}")

        per_contig = read_per_contig_file(args.per_contig) if args.per_contig else None
        supplied = supplied_seqkit or supplied_totals.get(args.sample)
        values, bases_total = assess_one(fields, supplied, per_contig)

        row = [args.sample, args.set_label] + [str(values[c]) for c in OUT_COLUMNS[2:]]
        text = "\t".join(OUT_COLUMNS) + "\n" + "\t".join(row) + "\n"
        if args.out:
            with open(args.out, "w") as fh:
                fh.write(text)
        else:
            sys.stdout.write(text)

        if args.per_genome_out:
            if not per_contig:
                sys.exit("--per-genome-out requires --per-contig")
            if not args.contig_map:
                sys.exit("--per-genome-out requires --contig-map")
            contig_genome = read_contig_map(args.contig_map)
            per_genome = aggregate_per_genome(per_contig, contig_genome)
            write_per_genome(args.per_genome_out, per_genome, args.sample,
                              args.set_label, bases_total)
        return

    # Standalone mode: one or more bare BAMs, upstream-compatible CLI.
    bams = collect_bams(args.bams, args.bam_dir)
    if not bams:
        sys.exit("no BAM files given (pass paths/globs/dirs or --bam-dir)")

    out_lines = ["\t".join(OUT_COLUMNS)]
    n_ok = 0
    for bam in bams:
        if not os.path.isfile(bam):
            warn(f"not a file, skipping: {bam}")
            continue
        sid = args.sample or sample_name(bam, args.strip_suffix)
        try:
            fields = run_samtools_stats(samtools, bam, args.threads)
            per_contig = compute_per_contig_from_bam(samtools, bam, args.threads)
        except RuntimeError as e:
            warn(f"skipping {bam}: {e}")
            continue
        if "reads_total" not in fields:
            warn(f"skipping {bam}: no summary numbers from samtools stats "
                 f"(empty or non-BAM?)")
            continue

        supplied = supplied_seqkit or supplied_totals.get(sid)
        values, _ = assess_one(fields, supplied, per_contig)
        row = [sid, args.set_label] + [str(values[c]) for c in OUT_COLUMNS[2:]]
        out_lines.append("\t".join(row))
        n_ok += 1

    if n_ok == 0:
        sys.exit("no BAMs could be processed")

    text = "\n".join(out_lines) + "\n"
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text)
        print(f"[assess_mapping] wrote {n_ok} row(s) to {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
