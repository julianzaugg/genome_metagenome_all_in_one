#!/usr/bin/env python
"""
Assemble a per-domain GTDB-Tk marker-protein MSA for tree building.

Re-implements, as one self-contained step, the logic of the example scripts
get_closest_reference_genomes.sh / get_related_reference_genomes.sh /
get_closest_leaves.py:

  1. Pick which user genomes to place: the supplied representative genomes,
     optionally filtered by CheckM2 completeness/contamination.
  2. Select GTDB reference genomes by:
       - closest-by-topology  (N closest reference leaves per genome in the
         GTDB-Tk classify tree; user genomes are never chosen as references)
       - related              (one reference per order that shares a bin's
         order but not its family)
       - explicit accessions  (a user-supplied list, added verbatim)
  3. Write <domain>.marker_msa.fasta = aligned user genomes (header carries the
     GTDB lineage) + aligned reference genomes, pulled from the GTDB-Tk MSAs.

Domains (bac120 / ar53) are handled independently; a domain with no user
genomes is skipped. Uses only the Python standard library plus ete3, so it has
no dependency on GNU awk/grep/zcat being present in the container.
"""

import argparse
import csv
import glob
import gzip
import os
import re
import sys

REF_PREFIXES = ("GB_", "RS_")


def eprint(*a):
    print(*a, file=sys.stderr)


def strip_ext(name):
    """Genome id from a filename / CheckM2 'Name' (drop a single extension)."""
    base = os.path.basename(name)
    return re.sub(r"\.(fa|fna|fasta)$", "", base)


def display_id(name, prefix):
    """Original genome name for display (drop the reserved reference prefix)."""
    return name[len(prefix):] if prefix and name.startswith(prefix) else name


def read_singleline_fasta(path, opener):
    """Yield (id, full_header, seq) from a FASTA whose records are single-line.

    GTDB-Tk MSAs are single-line per sequence; we still tolerate wrapped
    sequences by joining continuation lines.
    """
    hid = full = None
    seq = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hid is not None:
                    yield hid, full, "".join(seq)
                full = line[1:]
                hid = full.split()[0] if full.strip() else ""
                seq = []
            else:
                seq.append(line)
        if hid is not None:
            yield hid, full, "".join(seq)


def lineage_token(lineage, rank):
    """Return the e.g. 'f__Name' token from a GTDB lineage, '' if absent/blank."""
    for part in lineage.split(";"):
        part = part.strip()
        if part.startswith(rank) and len(part) > len(rank):
            return part
    return ""


def first_glob(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


def closest_reference_leaves(tree_file, query_ids, user_ids, n, exclude_pat):
    """For each query genome present in the tree, return its N closest leaves.

    User genomes are never returned (they cannot be references). Mirrors
    get_closest_leaves.py but loads the tree once and excludes the whole user
    set rather than relying on a single regex.
    """
    from ete3 import Tree

    tree = Tree(tree_file, format=1, quoted_node_names=True)
    pat = re.compile(exclude_pat) if exclude_pat else None
    pairs = []
    for qid in query_ids:
        nodes = tree.search_nodes(name=qid)
        if not nodes:
            continue
        q = nodes[0]
        dists = []
        for leaf in tree.iter_leaves():
            if leaf.name == qid or leaf.name in user_ids:
                continue
            if pat and pat.search(leaf.name):
                continue
            dists.append((leaf.name, q.get_distance(leaf)))
        dists.sort(key=lambda x: x[1])
        for name, _ in dists[:n]:
            pairs.append((qid, name))
    return pairs


def select_related(ref_records, bin_orders, bin_families, per_order):
    """One reference per order that shares a bin's order but not its family.

    ref_records: dict id -> lineage. Deterministic (sorted) so runs reproduce,
    unlike the original `shuf`-based script.
    """
    # one candidate per family
    by_family = {}
    for rid in sorted(ref_records):
        lin = ref_records[rid]
        fam = lineage_token(lin, "f__")
        order = lineage_token(lin, "o__")
        if not fam or fam in by_family:
            continue
        by_family[fam] = (rid, order, fam)

    chosen = []
    counts = {}
    for rid, order, fam in sorted(by_family.values()):
        if order in bin_orders and fam not in bin_families:
            if counts.get(order, 0) < per_order:
                counts[order] = counts.get(order, 0) + 1
                chosen.append(rid)
    return chosen


def parse_summary(path):
    """user_genome -> classification from a GTDB-Tk summary.tsv."""
    lineages = {}
    if not path or not os.path.exists(path):
        return lineages
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return lineages
        try:
            gi = header.index("user_genome")
            ci = header.index("classification")
        except ValueError:
            gi, ci = 0, 1
        for row in reader:
            if len(row) > max(gi, ci):
                lineages[row[gi]] = row[ci]
    return lineages


def quality_pass(checkm2, min_comp, max_cont):
    """Set of genome ids passing CheckM2 thresholds (None => no filtering)."""
    if not checkm2 or not os.path.exists(checkm2) or os.path.getsize(checkm2) == 0:
        return None
    keep = set()
    with open(checkm2) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                comp = float(row["Completeness"])
                cont = float(row["Contamination"])
            except (KeyError, ValueError, TypeError):
                continue
            if comp >= min_comp and cont <= max_cont:
                keep.add(strip_ext(row["Name"]))
    return keep


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtdbtk-dir", required=True, help="GTDB-Tk classify_wf output directory")
    ap.add_argument("--genomes-dir", required=True, help="Dir of user genome FASTAs to place")
    ap.add_argument("--checkm2", default=None, help="CheckM2 quality_report.tsv (optional)")
    ap.add_argument("--accessions", default=None, help="File of reference accessions, one per line (optional)")
    ap.add_argument("--domains", default="bac120,ar53")
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--closest-n", type=int, default=2)
    ap.add_argument("--related-per-order", type=int, default=1)
    ap.add_argument("--use-closest", action="store_true")
    ap.add_argument("--use-related", action="store_true")
    ap.add_argument("--exclude-pattern", default="", help="Extra regex to exclude from closest leaves")
    ap.add_argument("--reference-prefix", default="",
                    help="Reserved prefix marking external reference genomes; they are always "
                         "placed (bypass the quality filter) and shown with the prefix stripped")
    args = ap.parse_args()

    rep_ids = {strip_ext(f) for f in os.listdir(args.genomes_dir)}
    qual = quality_pass(args.checkm2, args.min_completeness, args.max_contamination)
    kept = rep_ids if qual is None else (rep_ids & qual)
    # External reference genomes are deliberate anchors — always place them, regardless
    # of the CheckM2 quality filter (they may have no row in the MAG CheckM2 report).
    ref_ids = {r for r in rep_ids if args.reference_prefix and r.startswith(args.reference_prefix)}
    kept = kept | ref_ids
    eprint(f"[marker-tree] {len(rep_ids)} genomes supplied ({len(ref_ids)} references), "
           f"{len(kept)} kept after quality filter")

    user_accessions = []
    if args.accessions and os.path.exists(args.accessions):
        with open(args.accessions) as fh:
            user_accessions = [ln.strip() for ln in fh if ln.strip()]

    any_ref_mode = args.use_closest or args.use_related or bool(user_accessions)
    total_refs = 0

    for domain in [d.strip() for d in args.domains.split(",") if d.strip()]:
        align = first_glob(os.path.join(args.gtdbtk_dir, "align", f"*.{domain}.user_msa.fasta.gz"))
        full = first_glob(os.path.join(args.gtdbtk_dir, "align", f"*.{domain}.msa.fasta.gz"))
        summary = first_glob(os.path.join(args.gtdbtk_dir, "classify", f"*.{domain}.summary.tsv"))
        treef = first_glob(os.path.join(args.gtdbtk_dir, "classify", f"*.{domain}.classify.tree"))

        if not align or not os.path.exists(align):
            continue
        user = {hid: (full_h, seq) for hid, full_h, seq in read_singleline_fasta(align, gzip.open)}
        domain_kept = sorted(k for k in kept if k in user)
        if not domain_kept:
            eprint(f"[marker-tree] {domain}: no kept genomes in this domain, skipping")
            continue

        lineages = parse_summary(summary)
        out_records = []  # (header_without_gt, seq)
        for bid in domain_kept:
            lin = lineages.get(bid, "")
            # look up the lineage with the (possibly prefixed) id, but label with the original
            out_records.append((f"{display_id(bid, args.reference_prefix)} {lin}".strip(), user[bid][1]))

        # Reference genomes available in the combined MSA (user + references)
        full_recs = {}
        ref_lineage = {}
        if full and os.path.exists(full):
            for hid, full_h, seq in read_singleline_fasta(full, gzip.open):
                full_recs[hid] = seq
                if hid.startswith(REF_PREFIXES):
                    parts = full_h.split(None, 1)
                    ref_lineage[hid] = parts[1] if len(parts) > 1 else ""

        refs = set()

        if args.use_closest and treef and os.path.exists(treef):
            pairs = closest_reference_leaves(
                treef, domain_kept, set(user.keys()), args.closest_n, args.exclude_pattern
            )
            with open(f"{domain}.closest_references.tsv", "w") as o:
                for q, leaf in pairs:
                    o.write(f"{display_id(q, args.reference_prefix)}\t{leaf}\n")
                    if leaf.startswith(REF_PREFIXES):
                        refs.add(leaf)

        if args.use_related and ref_lineage:
            bin_orders = {lineage_token(lineages.get(b, ""), "o__") for b in domain_kept}
            bin_orders.discard("")
            bin_families = {lineage_token(lineages.get(b, ""), "f__") for b in domain_kept}
            bin_families.discard("")
            refs.update(select_related(ref_lineage, bin_orders, bin_families, args.related_per_order))

        for acc in user_accessions:
            if acc in full_recs:
                refs.add(acc)

        ref_present = sorted(r for r in refs if r in full_recs and r not in user)
        total_refs += len(ref_present)
        if ref_present:
            with open(f"{domain}.reference_genomes.tsv", "w") as o:
                for r in ref_present:
                    o.write(f"{r}\t{ref_lineage.get(r, '')}\n")
                    out_records.append((f"{r} {ref_lineage.get(r, '')}".strip(), full_recs[r]))

        with open(f"{domain}.marker_msa.fasta", "w") as o:
            for header, seq in out_records:
                o.write(f">{header}\n{seq}\n")
        eprint(f"[marker-tree] {domain}: {len(domain_kept)} genomes + {len(ref_present)} references")

    if any_ref_mode and total_refs == 0:
        sys.exit(
            "ERROR: reference selection was requested but no reference sequences "
            "(>GB_/>RS_ records) were found in the GTDB-Tk MSA. Check that "
            "gtdbtk.<domain>.msa.fasta.gz contains reference genomes, or disable "
            "reference selection to build a genomes-only tree."
        )


if __name__ == "__main__":
    main()
