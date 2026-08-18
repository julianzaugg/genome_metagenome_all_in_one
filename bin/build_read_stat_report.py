#!/usr/bin/env python3
"""
Build a per-sample read-stat report (TSV) from the pipeline's existing outputs.

Sources (each directory is optional; whatever is present shapes the columns):
  --seqkit-dir       SeqKit stats per QC stage:   {id}.{stage}.seqkit_stats.tsv
  --scaffold-dir     CoverM reads-vs-assembly:    {id}_counts.tsv
  --scaffold-sr-dir  CoverM short reads-vs-assembly (hybrid nanopore isolate)
  --repmag-dir       CoverM reads-vs-all-dereplicated-rep-MAGs: {id}_abundances.tsv
  --hq-dir           HQ representative MAG FASTAs: {basename}.fasta
  --hq-repmag-dir    CoverM reads-vs-HQ-rep-MAGs-only:          {id}_abundances.tsv
  --hq-derep-repmag-dir  CoverM reads-vs-HQ-first-then-dereplicated MAGs: {id}_abundances.tsv
  --hq-ref-repmag-dir    CoverM reads-vs-(HQ MAGs + reference genomes) dereplicated set: {id}_abundances.tsv
  --ws-derep-repmag-dir  CoverM reads-vs-within-sample(or -group) dereplicated bins: {id}_abundances.tsv
  --ws-hq-repmag-dir     CoverM reads-vs-within-sample(or -group) HQ-first-then-dereplicated MAGs: {id}_abundances.tsv
  --assess-dir           MAPPING_ASSESS per-sample bases-mapped metrics:
                         {id}.{slot}.mapping_assessment.tsv (one row each)
  --assess-genome-dir    MAPPING_ASSESS per-genome bases-mapped metrics:
                         {id}.{slot}.mapping_assessment_per_genome.tsv
  --bases-metrics        which of A/B/C to include as Bases_mapped_* columns
                         (default: A,B,C -- see docs/output.md for definitions)

Rules:
  - All *read* percentages are vs Raw_count (raw input reads); all *base*
    percentages are vs the raw input bases (the value GBbp is derived from).
  - Paired reads are counted forward + reverse (SeqKit runs on both files,
    one row per file, so summing num_seqs/sum_len over rows gives the
    combined total).
  - Stage columns are labelled by the tool that produced them.
  - "Dereplicated_Bins" is the full set of dereplicated cluster representatives
    (any quality); "MAGs" is reserved for the high-quality subset.
  - HQ MAGs get two distinct counts, because they differ when dereplication
    leaves redundant near-identical lower-quality bins alongside an HQ rep:
    competitive mapping against the full rep set then splits reads across
    those siblings, undercounting the HQ rep specifically.
      * Reads_mapped_HQ_MAGs_count        — subset of --repmag-dir's full-set
        mapping, summed over rows matching an --hq-dir basename.
      * Reads_mapped_HQ_MAGs_direct_count — a separate --hq-repmag-dir
        mapping run against the HQ MAGs only, no competing siblings.
    The same subset-vs-direct split applies to Bases_mapped_{A,B,C}_HQ_MAGs,
    summed from --assess-genome-dir's Dereplicated_Bins per-genome rows over
    the same --hq-dir basenames.
  - Reads_mapped_HQ_Derep_MAGs_count is a third HQ mapping (--hq-derep-repmag-dir)
    against a differently constructed set: HQ MAGs extracted from the FULL
    pre-dereplication bin set and THEN dereplicated (vs the sets above, which
    dereplicate first and then keep the HQ representatives).
  - Reads_mapped_HQ_Ref_MAGs_count is a fourth HQ mapping (--hq-ref-repmag-dir)
    against the HQ MAGs dereplicated TOGETHER with external reference genomes.
  - Reads_mapped_PerSample_Derep_MAGs_count / Reads_mapped_PerSample_HQ_MAGs_count
    map each sample's reads to bins dereplicated WITHIN that sample (or group)
    only, not across all samples (--ws-derep-repmag-dir / --ws-hq-repmag-dir).
  - Bases_mapped_A/B/C are the base-level counterparts of every Reads_mapped_*
    column above (same reference set, same alignment set): A = full read
    length of every mapped read, B = aligned bases (CIGAR M/I/=/X) of primary
    alignments only, C = aligned bases of primary + supplementary alignments
    (the honest number for long reads that split across a contig/genome
    boundary). See docs/output.md and bin/assess_mapping.py.

Two report shapes:
  metagenome  Sample_ID, GBbp, Raw_count, <stage>_count/percent...,
              (Reads_mapped_<Label>_count/percent, Bases_mapped_A/B/C_<Label>_count/percent)
              for each of Scaffolds, Dereplicated_Bins(+HQ_MAGs), HQ_MAGs_direct,
              HQ_Derep_MAGs, HQ_Ref_MAGs, PerSample_Derep_MAGs, PerSample_HQ_MAGs
  isolate     Sample_ID, GBbp, Raw_count, <stage>_count/percent...,
              Covered_fraction, Mean_coverage, Read_count, Read_count_percent,
              Bases_mapped_A/B/C_count/percent
              (+ _SR variants when a short-read mapping is present)
"""
import argparse
import glob
import os
import sys

# Raw baseline stage names (illumina vs nanopore); both map to "Raw".
RAW_STAGES = ("raw", "raw_long")
# Non-raw QC stages in pipeline order, with their tool-name labels.
STAGE_ORDER = ["fastp", "porechop", "fastplong", "cleanifier"]
STAGE_LABEL = {
    "fastp": "Fastp",
    "porechop": "Porechop",
    "fastplong": "Fastplong",
    "cleanifier": "Cleanifier",
}

# assess_mapping.py metric letter -> its column name in the per-sample /
# per-genome mapping_assessment TSVs.
BASES_METRIC_COLUMN = {
    "A": "bases_mapped_readlen",
    "B": "bases_mapped_cigar",
    "C": "bases_mapped_cigar_all",
}
PER_GENOME_METRIC_COLUMN = {"A": "bases_A", "B": "bases_B", "C": "bases_C"}
ALL_BASES_METRICS = ["A", "B", "C"]


def read_tsv(path):
    """Return (header_list, list_of_row_lists) for a tab-separated file."""
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def col_index(header, suffix):
    """Index of the first column whose name equals or ends with `suffix`."""
    for i, name in enumerate(header):
        if name == suffix or name.endswith(" " + suffix):
            return i
    return None


# ---------------------------------------------------------------------------
# SeqKit parsing
# ---------------------------------------------------------------------------
def parse_seqkit_dir(directory):
    """
    Return {sample_id: {stage: {'reads': int, 'bases': int}}}.

    Filenames look like `{id}.{stage}.seqkit_stats.tsv`; the id may itself
    contain dots, so strip the fixed suffix then split the stage off the end.
    """
    out = {}
    suffix = ".seqkit_stats.tsv"
    for path in sorted(glob.glob(os.path.join(directory, "*" + suffix))):
        name = os.path.basename(path)[: -len(suffix)]
        if "." not in name:
            continue
        sample_id, stage = name.rsplit(".", 1)
        header, data = read_tsv(path)
        n_idx = col_index(header, "num_seqs")
        b_idx = col_index(header, "sum_len")
        reads = bases = 0
        for row in data:
            if n_idx is not None and n_idx < len(row):
                reads += int(float(row[n_idx]))
            if b_idx is not None and b_idx < len(row):
                bases += int(float(row[b_idx]))
        out.setdefault(sample_id, {})[stage] = {"reads": reads, "bases": bases}
    return out


# ---------------------------------------------------------------------------
# CoverM parsing (reads-mapped / Count column)
# ---------------------------------------------------------------------------
def parse_coverm_dir(directory, suffix, hq_names=None):
    """
    Parse CoverM TSVs in `directory` named `{id}{suffix}`.

    Returns {sample_id: {'count', 'covered_fraction', 'mean', 'hq_count'}}.
    `count` sums the Count column across genome rows (skipping any "unmapped"
    row). For single-row outputs (whole-assembly mapping) covered_fraction and
    mean are taken from that row. When `hq_names` is given, `hq_count` sums Count
    only over rows whose Genome basename is in that set.
    """
    out = {}
    for path in sorted(glob.glob(os.path.join(directory, "*" + suffix))):
        sample_id = os.path.basename(path)[: -len(suffix)]
        header, data = read_tsv(path)
        if not header:
            out[sample_id] = {}
            continue
        c_idx = col_index(header, "Count")
        cf_idx = col_index(header, "Covered Fraction")
        m_idx = col_index(header, "Mean")
        total = hq_total = 0.0
        cov = mean = None
        rows_used = 0
        for row in data:
            genome = row[0] if row else ""
            if genome.strip().lower() == "unmapped":
                continue
            count = 0.0
            if c_idx is not None and c_idx < len(row):
                count = float(row[c_idx] or 0)
            total += count
            if hq_names is not None and strip_ext(genome) in hq_names:
                hq_total += count
            if rows_used == 0:
                if cf_idx is not None and cf_idx < len(row):
                    cov = float(row[cf_idx] or 0)
                if m_idx is not None and m_idx < len(row):
                    mean = float(row[m_idx] or 0)
            rows_used += 1
        rec = {"count": int(round(total))}
        if cov is not None:
            rec["covered_fraction"] = cov
        if mean is not None:
            rec["mean"] = mean
        if hq_names is not None:
            rec["hq_count"] = int(round(hq_total))
        out[sample_id] = rec
    return out


def strip_ext(name):
    base = os.path.basename(name)
    for ext in (".fasta", ".fa", ".fna"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def hq_basenames(directory):
    names = set()
    for path in glob.glob(os.path.join(directory, "*")):
        if os.path.isfile(path) and not path.endswith("README.txt"):
            names.add(strip_ext(path))
    return names


# ---------------------------------------------------------------------------
# MAPPING_ASSESS parsing (bases-mapped / percent-of-sequenced-bases)
# ---------------------------------------------------------------------------
def parse_assess_dir(directory):
    """
    Parse `{id}.{slot}.mapping_assessment.tsv` files (one data row each).

    Returns {sample_id: {slot: {metric_letter: bases_count}}}.
    """
    out = {}
    suffix = ".mapping_assessment.tsv"
    for path in sorted(glob.glob(os.path.join(directory, "*" + suffix))):
        name = os.path.basename(path)[: -len(suffix)]
        if "." not in name:
            continue
        sample_id, slot = name.rsplit(".", 1)
        header, data = read_tsv(path)
        if not header or not data:
            continue
        row = data[0]
        rec = {}
        for letter, col in BASES_METRIC_COLUMN.items():
            idx = col_index(header, col)
            if idx is not None and idx < len(row) and row[idx] != "":
                rec[letter] = int(float(row[idx]))
        out.setdefault(sample_id, {})[slot] = rec
    return out


def parse_assess_genome_dir(directory):
    """
    Parse `{id}.{slot}.mapping_assessment_per_genome.tsv` files.

    Returns {sample_id: {slot: {genome: {metric_letter: bases_count}}}}.
    """
    out = {}
    suffix = ".mapping_assessment_per_genome.tsv"
    for path in sorted(glob.glob(os.path.join(directory, "*" + suffix))):
        name = os.path.basename(path)[: -len(suffix)]
        if "." not in name:
            continue
        sample_id, slot = name.rsplit(".", 1)
        header, data = read_tsv(path)
        if not header:
            continue
        g_idx = col_index(header, "genome")
        metric_idx = {letter: col_index(header, col)
                      for letter, col in PER_GENOME_METRIC_COLUMN.items()}
        genomes = {}
        for row in data:
            if g_idx is None or g_idx >= len(row):
                continue
            genome = row[g_idx]
            rec = {}
            for letter, idx in metric_idx.items():
                if idx is not None and idx < len(row) and row[idx] != "":
                    rec[letter] = int(float(row[idx]))
            genomes[genome] = rec
        out.setdefault(sample_id, {})[slot] = genomes
    return out


def sum_bases_over_genomes(genome_table, hq_names, letter):
    """Sum one metric over genome rows whose (extension-stripped) name is in
    hq_names -- the base-level counterpart of parse_coverm_dir's hq_count."""
    total = 0
    for genome, rec in genome_table.items():
        if strip_ext(genome) in hq_names and letter in rec:
            total += rec[letter]
    return total


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------
def pct(count, raw):
    if count is None or not raw:
        return ""
    return f"{count / raw * 100:.2f}"


def raw_record(stages):
    """Return the raw stage record for a sample, or None."""
    for s in RAW_STAGES:
        if s in stages:
            return stages[s]
    return None


def bases_cols(header_or_row, label, assess_slot, metrics, is_header, sample_bases=None,
               hq_names=None, genome_table=None):
    """
    Build the Bases_mapped_{letter}_{label}_count/_percent header cells, or
    (for a data row) their values for one sample.

    assess_slot: {letter: bases_count} for the sample+slot (from --assess-dir),
    or None if that slot's assessment is missing for this sample.
    When hq_names + genome_table are given, values are instead summed from the
    per-genome table over genomes in hq_names (the HQ-subset column).
    """
    cells = []
    for letter in metrics:
        if is_header:
            cells += [f"Bases_mapped_{letter}_{label}_count", f"Bases_mapped_{letter}_{label}_percent"]
        else:
            if hq_names is not None:
                count = sum_bases_over_genomes(genome_table or {}, hq_names, letter) if genome_table else None
            else:
                count = (assess_slot or {}).get(letter)
            cells += ["" if count is None else str(count), pct(count, sample_bases)]
    return cells


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------
def build(mode, seqkit, scaffold, scaffold_sr, repmag, hq_names, has_hq,
          repmag_direct, repmag_hq_derep, repmag_hq_ref,
          repmag_ws_derep, repmag_ws_hq, assess, assess_genome, metrics):
    samples = sorted(seqkit)

    present_stages = set()
    for stages in seqkit.values():
        present_stages.update(stages)
    qc_stages = [s for s in STAGE_ORDER if s in present_stages]

    header = ["Sample_ID", "GBbp", "Raw_count"]
    for s in qc_stages:
        header += [f"{STAGE_LABEL[s]}_count", f"{STAGE_LABEL[s]}_percent"]

    # (slot label used for --assess-dir lookups, report column label, source dict)
    if mode == "metagenome":
        slots = [
            ("Scaffolds", "Scaffolds", scaffold),
            ("Dereplicated_Bins", "Dereplicated_Bins", repmag),
            ("HQ_MAGs_direct", "HQ_MAGs_direct", repmag_direct),
            ("HQ_Derep_MAGs", "HQ_Derep_MAGs", repmag_hq_derep),
            ("HQ_Ref_MAGs", "HQ_Ref_MAGs", repmag_hq_ref),
            ("PerSample_Derep_MAGs", "PerSample_Derep_MAGs", repmag_ws_derep),
            ("PerSample_HQ_MAGs", "PerSample_HQ_MAGs", repmag_ws_hq),
        ]
        if scaffold:
            header += ["Reads_mapped_Scaffolds_count", "Reads_mapped_Scaffolds_percent"]
            header += bases_cols(None, "Scaffolds", None, metrics, True)
        if repmag:
            header += ["Reads_mapped_Dereplicated_Bins_count", "Reads_mapped_Dereplicated_Bins_percent"]
            header += bases_cols(None, "Dereplicated_Bins", None, metrics, True)
            if has_hq:
                header += ["Reads_mapped_HQ_MAGs_count", "Reads_mapped_HQ_MAGs_percent"]
                header += bases_cols(None, "HQ_MAGs", None, metrics, True)
        if repmag_direct:
            header += ["Reads_mapped_HQ_MAGs_direct_count", "Reads_mapped_HQ_MAGs_direct_percent"]
            header += bases_cols(None, "HQ_MAGs_direct", None, metrics, True)
        if repmag_hq_derep:
            header += ["Reads_mapped_HQ_Derep_MAGs_count", "Reads_mapped_HQ_Derep_MAGs_percent"]
            header += bases_cols(None, "HQ_Derep_MAGs", None, metrics, True)
        if repmag_hq_ref:
            header += ["Reads_mapped_HQ_Ref_MAGs_count", "Reads_mapped_HQ_Ref_MAGs_percent"]
            header += bases_cols(None, "HQ_Ref_MAGs", None, metrics, True)
        if repmag_ws_derep:
            header += ["Reads_mapped_PerSample_Derep_MAGs_count", "Reads_mapped_PerSample_Derep_MAGs_percent"]
            header += bases_cols(None, "PerSample_Derep_MAGs", None, metrics, True)
        if repmag_ws_hq:
            header += ["Reads_mapped_PerSample_HQ_MAGs_count", "Reads_mapped_PerSample_HQ_MAGs_percent"]
            header += bases_cols(None, "PerSample_HQ_MAGs", None, metrics, True)
    else:  # isolate
        if scaffold:
            header += ["Covered_fraction", "Mean_coverage", "Read_count", "Read_count_percent"]
            header += bases_cols(None, "Scaffolds", None, metrics, True)
        if scaffold_sr:
            header += ["Covered_fraction_SR", "Mean_coverage_SR", "Read_count_SR", "Read_count_SR_percent"]
            header += bases_cols(None, "Scaffolds_SR", None, metrics, True)

    lines = ["\t".join(header)]
    for sid in samples:
        stages = seqkit[sid]
        raw = raw_record(stages)
        raw_reads = raw["reads"] if raw else 0
        raw_bases = raw["bases"] if raw else 0
        gbbp = f"{raw_bases / 1e9:.5f}" if raw else ""
        row = [sid, gbbp, str(raw_reads) if raw else ""]

        for s in qc_stages:
            if s in stages:
                c = stages[s]["reads"]
                row += [str(c), pct(c, raw_reads)]
            else:
                row += ["", ""]

        sample_assess = assess.get(sid, {})
        sample_assess_genome = assess_genome.get(sid, {})

        if mode == "metagenome":
            if scaffold:
                rec = scaffold.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                row += bases_cols(None, "Scaffolds", sample_assess.get("Scaffolds"), metrics, False, raw_bases)
            if repmag:
                rec = repmag.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                row += bases_cols(None, "Dereplicated_Bins", sample_assess.get("Dereplicated_Bins"), metrics, False, raw_bases)
                if has_hq:
                    h = rec.get("hq_count")
                    row += ["" if h is None else str(h), pct(h, raw_reads) if h is not None else ""]
                    row += bases_cols(None, "HQ_MAGs", None, metrics, False, raw_bases,
                                       hq_names=hq_names,
                                       genome_table=sample_assess_genome.get("Dereplicated_Bins"))
            if repmag_direct:
                rec = repmag_direct.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                row += bases_cols(None, "HQ_MAGs_direct", sample_assess.get("HQ_MAGs_direct"), metrics, False, raw_bases)
            if repmag_hq_derep:
                rec = repmag_hq_derep.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                row += bases_cols(None, "HQ_Derep_MAGs", sample_assess.get("HQ_Derep_MAGs"), metrics, False, raw_bases)
            if repmag_hq_ref:
                rec = repmag_hq_ref.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                row += bases_cols(None, "HQ_Ref_MAGs", sample_assess.get("HQ_Ref_MAGs"), metrics, False, raw_bases)
            if repmag_ws_derep:
                rec = repmag_ws_derep.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                row += bases_cols(None, "PerSample_Derep_MAGs", sample_assess.get("PerSample_Derep_MAGs"), metrics, False, raw_bases)
            if repmag_ws_hq:
                rec = repmag_ws_hq.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                row += bases_cols(None, "PerSample_HQ_MAGs", sample_assess.get("PerSample_HQ_MAGs"), metrics, False, raw_bases)
        else:
            if scaffold:
                rec = scaffold.get(sid, {})
                row += _isolate_cols(rec, raw_reads)
                row += bases_cols(None, "Scaffolds", sample_assess.get("Scaffolds"), metrics, False, raw_bases)
            if scaffold_sr:
                rec = scaffold_sr.get(sid, {})
                row += _isolate_cols(rec, raw_reads)
                row += bases_cols(None, "Scaffolds_SR", sample_assess.get("Scaffolds_SR"), metrics, False, raw_bases)

        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def _isolate_cols(rec, raw_reads):
    cov = rec.get("covered_fraction")
    mean = rec.get("mean")
    cnt = rec.get("count")
    return [
        "" if cov is None else f"{cov:.4f}",
        "" if mean is None else f"{mean:.4f}",
        "" if cnt is None else str(cnt),
        pct(cnt, raw_reads) if cnt is not None else "",
    ]


def maybe_dir(path):
    return path if path and os.path.isdir(path) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["metagenome", "isolate"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--seqkit-dir", required=True)
    ap.add_argument("--scaffold-dir", default="")
    ap.add_argument("--scaffold-sr-dir", default="")
    ap.add_argument("--repmag-dir", default="")
    ap.add_argument("--hq-dir", default="")
    ap.add_argument("--hq-repmag-dir", default="")
    ap.add_argument("--hq-derep-repmag-dir", default="")
    ap.add_argument("--hq-ref-repmag-dir", default="")
    ap.add_argument("--ws-derep-repmag-dir", default="")
    ap.add_argument("--ws-hq-repmag-dir", default="")
    ap.add_argument("--assess-dir", default="")
    ap.add_argument("--assess-genome-dir", default="")
    ap.add_argument("--bases-metrics", default="A,B,C",
                    help="comma-separated subset of A,B,C to include as "
                         "Bases_mapped_* columns (default: A,B,C)")
    args = ap.parse_args()

    if not os.path.isdir(args.seqkit_dir):
        sys.exit(f"seqkit dir not found: {args.seqkit_dir}")
    seqkit = parse_seqkit_dir(args.seqkit_dir)

    metrics = [m.strip().upper() for m in args.bases_metrics.split(",") if m.strip()]
    bad = [m for m in metrics if m not in ALL_BASES_METRICS]
    if bad:
        sys.exit(f"--bases-metrics: unknown metric(s) {bad}; choose from A, B, C")

    hq_dir = maybe_dir(args.hq_dir)
    hq_names = hq_basenames(hq_dir) if hq_dir else set()
    has_hq = hq_dir is not None and bool(hq_names)

    scaffold_dir = maybe_dir(args.scaffold_dir)
    scaffold = parse_coverm_dir(scaffold_dir, "_counts.tsv") if scaffold_dir else {}
    sr_dir = maybe_dir(args.scaffold_sr_dir)
    scaffold_sr = parse_coverm_dir(sr_dir, "_counts.tsv") if sr_dir else {}
    repmag_dir = maybe_dir(args.repmag_dir)
    repmag = parse_coverm_dir(repmag_dir, "_abundances.tsv",
                              hq_names=hq_names if has_hq else None) if repmag_dir else {}

    hq_repmag_dir = maybe_dir(args.hq_repmag_dir)
    repmag_direct = parse_coverm_dir(hq_repmag_dir, "_abundances.tsv") if hq_repmag_dir else {}

    hq_derep_dir = maybe_dir(args.hq_derep_repmag_dir)
    repmag_hq_derep = parse_coverm_dir(hq_derep_dir, "_abundances.tsv") if hq_derep_dir else {}

    hq_ref_dir = maybe_dir(args.hq_ref_repmag_dir)
    repmag_hq_ref = parse_coverm_dir(hq_ref_dir, "_abundances.tsv") if hq_ref_dir else {}

    ws_derep_dir = maybe_dir(args.ws_derep_repmag_dir)
    repmag_ws_derep = parse_coverm_dir(ws_derep_dir, "_abundances.tsv") if ws_derep_dir else {}

    ws_hq_dir = maybe_dir(args.ws_hq_repmag_dir)
    repmag_ws_hq = parse_coverm_dir(ws_hq_dir, "_abundances.tsv") if ws_hq_dir else {}

    assess_dir = maybe_dir(args.assess_dir)
    assess = parse_assess_dir(assess_dir) if assess_dir else {}
    assess_genome_dir = maybe_dir(args.assess_genome_dir)
    assess_genome = parse_assess_genome_dir(assess_genome_dir) if assess_genome_dir else {}

    report = build(args.mode, seqkit, scaffold, scaffold_sr,
                   repmag, hq_names, has_hq, repmag_direct, repmag_hq_derep, repmag_hq_ref,
                   repmag_ws_derep, repmag_ws_hq, assess, assess_genome, metrics)
    with open(args.out, "w") as fh:
        fh.write(report)


if __name__ == "__main__":
    main()
