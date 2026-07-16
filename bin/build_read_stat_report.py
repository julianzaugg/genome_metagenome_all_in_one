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

Rules:
  - All percentages are vs Raw (raw input reads). Raw has no percent.
  - Paired reads are counted forward + reverse (SeqKit runs on both files,
    one row per file, so summing num_seqs over rows gives the combined total).
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
  - Reads_mapped_HQ_Derep_MAGs_count is a third HQ mapping (--hq-derep-repmag-dir)
    against a differently constructed set: HQ MAGs extracted from the FULL
    pre-dereplication bin set and THEN dereplicated (vs the sets above, which
    dereplicate first and then keep the HQ representatives).
  - Reads_mapped_HQ_Ref_MAGs_count is a fourth HQ mapping (--hq-ref-repmag-dir)
    against the HQ MAGs dereplicated TOGETHER with external reference genomes.
  - Reads_mapped_PerSample_Derep_MAGs_count / Reads_mapped_PerSample_HQ_MAGs_count
    map each sample's reads to bins dereplicated WITHIN that sample (or group)
    only, not across all samples (--ws-derep-repmag-dir / --ws-hq-repmag-dir).

Two report shapes:
  metagenome  Sample_ID, GBbp, Raw_count, <stage>_count/percent...,
              Reads_mapped_{Scaffolds,Dereplicated_Bins,HQ_MAGs,HQ_MAGs_direct,HQ_Derep_MAGs}_count/percent
  isolate     Sample_ID, GBbp, Raw_count, <stage>_count/percent...,
              Covered_fraction, Mean_coverage, Read_count, Read_count_percent
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
# CoverM parsing
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
# Formatting helpers
# ---------------------------------------------------------------------------
def pct(count, raw):
    if not raw:
        return ""
    return f"{count / raw * 100:.2f}"


def raw_record(stages):
    """Return the raw stage record for a sample, or None."""
    for s in RAW_STAGES:
        if s in stages:
            return stages[s]
    return None


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------
def build(mode, seqkit, scaffold, scaffold_sr, repmag, hq_names, has_hq,
          repmag_direct, repmag_hq_derep, repmag_hq_ref,
          repmag_ws_derep, repmag_ws_hq):
    samples = sorted(seqkit)

    present_stages = set()
    for stages in seqkit.values():
        present_stages.update(stages)
    qc_stages = [s for s in STAGE_ORDER if s in present_stages]

    header = ["Sample_ID", "GBbp", "Raw_count"]
    for s in qc_stages:
        header += [f"{STAGE_LABEL[s]}_count", f"{STAGE_LABEL[s]}_percent"]

    if mode == "metagenome":
        if scaffold:
            header += ["Reads_mapped_Scaffolds_count", "Reads_mapped_Scaffolds_percent"]
        if repmag:
            header += ["Reads_mapped_Dereplicated_Bins_count", "Reads_mapped_Dereplicated_Bins_percent"]
            if has_hq:
                header += ["Reads_mapped_HQ_MAGs_count", "Reads_mapped_HQ_MAGs_percent"]
        if repmag_direct:
            header += ["Reads_mapped_HQ_MAGs_direct_count", "Reads_mapped_HQ_MAGs_direct_percent"]
        if repmag_hq_derep:
            header += ["Reads_mapped_HQ_Derep_MAGs_count", "Reads_mapped_HQ_Derep_MAGs_percent"]
        if repmag_hq_ref:
            header += ["Reads_mapped_HQ_Ref_MAGs_count", "Reads_mapped_HQ_Ref_MAGs_percent"]
        if repmag_ws_derep:
            header += ["Reads_mapped_PerSample_Derep_MAGs_count", "Reads_mapped_PerSample_Derep_MAGs_percent"]
        if repmag_ws_hq:
            header += ["Reads_mapped_PerSample_HQ_MAGs_count", "Reads_mapped_PerSample_HQ_MAGs_percent"]
    else:  # isolate
        if scaffold:
            header += ["Covered_fraction", "Mean_coverage", "Read_count", "Read_count_percent"]
        if scaffold_sr:
            header += ["Covered_fraction_SR", "Mean_coverage_SR", "Read_count_SR", "Read_count_SR_percent"]

    lines = ["\t".join(header)]
    for sid in samples:
        stages = seqkit[sid]
        raw = raw_record(stages)
        raw_reads = raw["reads"] if raw else 0
        gbbp = f"{raw['bases'] / 1e9:.5f}" if raw else ""
        row = [sid, gbbp, str(raw_reads) if raw else ""]

        for s in qc_stages:
            if s in stages:
                c = stages[s]["reads"]
                row += [str(c), pct(c, raw_reads)]
            else:
                row += ["", ""]

        if mode == "metagenome":
            if scaffold:
                rec = scaffold.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
            if repmag:
                rec = repmag.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
                if has_hq:
                    h = rec.get("hq_count")
                    row += ["" if h is None else str(h), pct(h, raw_reads) if h is not None else ""]
            if repmag_direct:
                rec = repmag_direct.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
            if repmag_hq_derep:
                rec = repmag_hq_derep.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
            if repmag_hq_ref:
                rec = repmag_hq_ref.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
            if repmag_ws_derep:
                rec = repmag_ws_derep.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
            if repmag_ws_hq:
                rec = repmag_ws_hq.get(sid, {})
                c = rec.get("count")
                row += ["" if c is None else str(c), pct(c, raw_reads) if c is not None else ""]
        else:
            if scaffold:
                rec = scaffold.get(sid, {})
                row += _isolate_cols(rec, raw_reads)
            if scaffold_sr:
                rec = scaffold_sr.get(sid, {})
                row += _isolate_cols(rec, raw_reads)

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
    ap = argparse.ArgumentParser(description=__doc__)
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
    args = ap.parse_args()

    if not os.path.isdir(args.seqkit_dir):
        sys.exit(f"seqkit dir not found: {args.seqkit_dir}")
    seqkit = parse_seqkit_dir(args.seqkit_dir)

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

    report = build(args.mode, seqkit, scaffold, scaffold_sr,
                   repmag, hq_names, has_hq, repmag_direct, repmag_hq_derep, repmag_hq_ref,
                   repmag_ws_derep, repmag_ws_hq)
    with open(args.out, "w") as fh:
        fh.write(report)


if __name__ == "__main__":
    main()
