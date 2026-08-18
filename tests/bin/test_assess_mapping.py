#!/usr/bin/env python3

import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "assess_mapping.py"
BAM_DIR = REPO / "tests" / "data" / "bams"

spec = importlib.util.spec_from_file_location("assess_mapping", SCRIPT)
assess_mapping = importlib.util.module_from_spec(spec)
sys.modules["assess_mapping"] = assess_mapping
spec.loader.exec_module(assess_mapping)

HAVE_SAMTOOLS = shutil.which("samtools") is not None


def read_tsv(path):
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return rows[0], rows[1:]


class ParseSnTest(unittest.TestCase):
    def test_parses_known_fields_and_ignores_others(self):
        text = (
            "SN\traw total sequences:\t3\t# excluding supplementary and secondary reads\n"
            "SN\treads mapped:\t2\n"
            "SN\treads unmapped:\t1\n"
            "SN\tbases mapped:\t2500\t# ignores clipping\n"
            "SN\tbases mapped (cigar):\t2300\t# more accurate\n"
            "SN\ttotal length:\t7500\n"
            "SN\terror rate:\t0.0\t# not tracked\n"
        )
        fields = assess_mapping.parse_sn(text)
        self.assertEqual(fields, {
            "reads_total": 3, "reads_mapped": 2, "reads_unmapped": 1,
            "bases_mapped_readlen": 2500, "bases_mapped_cigar": 2300,
            "bases_total": 7500,
        })

    def test_ignores_non_sn_lines(self):
        text = "CHK\tsome\tother\tline\nSN\treads mapped:\t5\n"
        self.assertEqual(assess_mapping.parse_sn(text), {"reads_mapped": 5})


class CigarAlignedBasesTest(unittest.TestCase):
    def test_simple_match(self):
        self.assertEqual(assess_mapping.cigar_aligned_bases("500M"), 500)

    def test_soft_clip_excluded(self):
        self.assertEqual(assess_mapping.cigar_aligned_bases("1800M200S"), 1800)

    def test_hard_clip_excluded_insertion_included(self):
        # 15000M15000H -> 15000 (H carries no read bases at all, M counts)
        self.assertEqual(assess_mapping.cigar_aligned_bases("15000M15000H"), 15000)
        # complementary supplementary record for the same split read
        self.assertEqual(assess_mapping.cigar_aligned_bases("15000H15000M"), 15000)

    def test_deletion_excluded(self):
        # D consumes reference, not read bases -> not counted
        self.assertEqual(assess_mapping.cigar_aligned_bases("50M10D50M"), 100)

    def test_star_cigar_is_zero(self):
        self.assertEqual(assess_mapping.cigar_aligned_bases("*"), 0)


class ResolveDenominatorsTest(unittest.TestCase):
    def test_bam_source_when_unmapped_present(self):
        fields = {"reads_unmapped": 1, "reads_total": 3, "bases_total": 7500}
        reads_total, bases_total, source, note = assess_mapping.resolve_denominators(fields, None)
        self.assertEqual((reads_total, bases_total, source, note), (3, 7500, "bam", ""))

    def test_supplied_source_when_no_unmapped(self):
        fields = {"reads_unmapped": 0, "reads_total": 2, "bases_total": 4000}
        supplied = {"total_reads": 4, "total_bases": 10000}
        reads_total, bases_total, source, note = assess_mapping.resolve_denominators(fields, supplied)
        self.assertEqual((reads_total, bases_total, source), (4, 10000, "supplied"))

    def test_bam_mapped_only_warns_when_nothing_supplied(self):
        fields = {"reads_unmapped": 0, "reads_total": 2, "bases_total": 4000}
        reads_total, bases_total, source, note = assess_mapping.resolve_denominators(fields, None)
        self.assertEqual((reads_total, bases_total, source), (2, 4000, "bam(mapped-only)"))
        self.assertTrue(note)


class AssessOneTest(unittest.TestCase):
    def test_metric_c_falls_back_to_b_without_per_contig(self):
        fields = {"reads_mapped": 2, "bases_mapped_readlen": 2500,
                  "bases_mapped_cigar": 2300, "reads_total": 3,
                  "reads_unmapped": 1, "bases_total": 7500}
        values, bases_total = assess_mapping.assess_one(fields, None, per_contig=None)
        self.assertEqual(values["bases_mapped_cigar_all"], 2300)
        self.assertEqual(values["reads_supplementary"], 0)
        self.assertEqual(bases_total, 7500)

    def test_metric_c_sums_per_contig_including_supplementary(self):
        fields = {"reads_mapped": 2, "bases_mapped_readlen": 9500,
                  "bases_mapped_cigar": 5503, "reads_total": 3,
                  "reads_unmapped": 0, "bases_total": 11500}
        per_contig = {
            "genomeA": {"reads_primary": 2, "reads_supplementary": 0,
                        "bases_A": 9500, "bases_B": 5503, "bases_C": 5503},
            "genomeB": {"reads_primary": 0, "reads_supplementary": 1,
                        "bases_A": 0, "bases_B": 0, "bases_C": 4000},
        }
        supplied = {"total_reads": 3, "total_bases": 11500}
        values, bases_total = assess_mapping.assess_one(fields, supplied, per_contig)
        self.assertEqual(values["bases_mapped_cigar_all"], 9503)
        self.assertEqual(values["reads_supplementary"], 1)
        self.assertEqual(values["pct_bases_mapped_C"], "82.63")
        self.assertEqual(values["total_source"], "supplied")

    def test_self_check_warns_on_mismatch(self):
        fields = {"reads_mapped": 1, "bases_mapped_readlen": 100,
                  "bases_mapped_cigar": 90, "reads_total": 1,
                  "reads_unmapped": 0, "bases_total": 100}
        # deliberately wrong per-contig sum (bases_A should be 100, not 50)
        per_contig = {"c1": {"reads_primary": 1, "reads_supplementary": 0,
                              "bases_A": 50, "bases_B": 90, "bases_C": 90}}
        buf = []
        orig_warn = assess_mapping.warn
        assess_mapping.warn = lambda msg: buf.append(msg)
        try:
            assess_mapping.assess_one(fields, {"total_reads": 1, "total_bases": 100}, per_contig)
        finally:
            assess_mapping.warn = orig_warn
        self.assertTrue(any("self-check failed" in m for m in buf))


class AggregatePerGenomeTest(unittest.TestCase):
    def test_groups_by_contig_map_and_flags_unassigned(self):
        per_contig = {
            "mag1~contigA": {"reads_primary": 5, "reads_supplementary": 1,
                              "bases_A": 500, "bases_B": 450, "bases_C": 480},
            "mag1~contigB": {"reads_primary": 2, "reads_supplementary": 0,
                              "bases_A": 200, "bases_B": 180, "bases_C": 180},
            "unknown~contigX": {"reads_primary": 1, "reads_supplementary": 0,
                                 "bases_A": 100, "bases_B": 90, "bases_C": 90},
        }
        contig_map = {"mag1~contigA": "mag1", "mag1~contigB": "mag1"}
        buf = []
        orig_warn = assess_mapping.warn
        assess_mapping.warn = lambda msg: buf.append(msg)
        try:
            per_genome = assess_mapping.aggregate_per_genome(per_contig, contig_map)
        finally:
            assess_mapping.warn = orig_warn
        self.assertEqual(per_genome["mag1"]["reads_primary"], 7)
        self.assertEqual(per_genome["mag1"]["bases_A"], 700)
        self.assertEqual(per_genome["unassigned"]["reads_primary"], 1)
        self.assertTrue(any("unassigned" in m for m in buf))


@unittest.skipUnless(HAVE_SAMTOOLS, "samtools not on PATH")
class StandaloneModeRegressionTest(unittest.TestCase):
    """Regression test against the upstream assess_long_read_mapping example
    fixtures -- metrics A/B and the read-count columns must match its
    documented expected_output.tsv exactly."""

    def test_sampleA_totals_from_bam(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.tsv"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(BAM_DIR / "sampleA.bam"), "-o", str(out)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            header, rows = read_tsv(out)
            row = dict(zip(header, rows[0]))
            self.assertEqual(row["reads_mapped"], "2")
            self.assertEqual(row["reads_total"], "3")
            self.assertEqual(row["pct_reads_mapped"], "66.67")
            self.assertEqual(row["bases_mapped_readlen"], "2500")
            self.assertEqual(row["pct_bases_mapped_A"], "33.33")
            self.assertEqual(row["bases_mapped_cigar"], "2300")
            self.assertEqual(row["pct_bases_mapped_B"], "30.67")
            self.assertEqual(row["bases_total"], "7500")
            self.assertEqual(row["total_source"], "bam")

    def test_sampleB_totals_from_supplied_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.tsv"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(BAM_DIR / "sampleB.bam"),
                 "--totals", str(BAM_DIR / "totals.tsv"), "-o", str(out)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            header, rows = read_tsv(out)
            row = dict(zip(header, rows[0]))
            self.assertEqual(row["reads_mapped"], "2")
            self.assertEqual(row["reads_total"], "4")
            self.assertEqual(row["pct_reads_mapped"], "50.00")
            self.assertEqual(row["bases_mapped_readlen"], "4000")
            self.assertEqual(row["pct_bases_mapped_A"], "40.00")
            self.assertEqual(row["bases_total"], "10000")
            self.assertEqual(row["total_source"], "supplied")

    def test_no_self_check_warning_on_real_bams(self):
        # sampleA/B have no supplementary alignments -> the pure-Python
        # per-contig fallback must agree with samtools stats exactly.
        for bam in ("sampleA.bam", "sampleB.bam"):
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(BAM_DIR / bam),
                 "--totals", str(BAM_DIR / "totals.tsv")],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("self-check failed", result.stderr)


@unittest.skipUnless(HAVE_SAMTOOLS, "samtools not on PATH")
class PipelineModeTest(unittest.TestCase):
    """Exercises --stats-file/--per-contig/--contig-map/--totals-from-seqkit,
    the interface MAPPING_ASSESS actually calls."""

    def test_pipeline_mode_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            stats_file = tmp / "stats.txt"
            stats_file.write_text(
                "SN\traw total sequences:\t3\n"
                "SN\treads mapped:\t2\n"
                "SN\treads unmapped:\t0\n"
                "SN\tbases mapped:\t9500\n"
                "SN\tbases mapped (cigar):\t5503\n"
                "SN\ttotal length:\t11500\n"
            )
            per_contig = tmp / "per_contig.tsv"
            per_contig.write_text(
                "contig\treads_primary\treads_supplementary\tbases_A\tbases_B\tbases_C\n"
                "genomeA~genomeA\t2\t0\t9500\t5503\t5503\n"
                "genomeB~genomeB\t0\t1\t0\t0\t4000\n"
            )
            contig_map = tmp / "contig_map.tsv"
            contig_map.write_text("genomeA~genomeA\tgenomeA\ngenomeB~genomeB\tgenomeB\n")
            seqkit = tmp / "seqkit.tsv"
            seqkit.write_text(
                "file\tformat\ttype\tnum_seqs\tsum_len\n"
                "reads.fastq\tFASTQ\tDNA\t3\t11500\n"
            )
            sample_out = tmp / "sample.tsv"
            genome_out = tmp / "genome.tsv"

            result = subprocess.run(
                [sys.executable, str(SCRIPT),
                 "--sample", "JZ_01", "--set", "HQ_Derep_MAGs",
                 "--stats-file", str(stats_file),
                 "--totals-from-seqkit", str(seqkit),
                 "--per-contig", str(per_contig),
                 "--contig-map", str(contig_map),
                 "--per-genome-out", str(genome_out),
                 "-o", str(sample_out)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("self-check failed", result.stderr)

            header, rows = read_tsv(sample_out)
            row = dict(zip(header, rows[0]))
            self.assertEqual(row["sample"], "JZ_01")
            self.assertEqual(row["set"], "HQ_Derep_MAGs")
            self.assertEqual(row["reads_mapped"], "2")
            self.assertEqual(row["reads_supplementary"], "1")
            self.assertEqual(row["bases_mapped_readlen"], "9500")
            self.assertEqual(row["bases_mapped_cigar"], "5503")
            self.assertEqual(row["bases_mapped_cigar_all"], "9503")
            self.assertEqual(row["pct_bases_mapped_C"], "82.63")
            self.assertEqual(row["total_source"], "supplied")

            g_header, g_rows = read_tsv(genome_out)
            genomes = {r[g_header.index("genome")]: r for r in g_rows}
            self.assertIn("genomeA", genomes)
            self.assertIn("genomeB", genomes)
            row_b = dict(zip(g_header, genomes["genomeB"]))
            self.assertEqual(row_b["bases_C"], "4000")
            self.assertEqual(row_b["reads_primary"], "0")
            self.assertEqual(row_b["reads_supplementary"], "1")

    def test_missing_stats_file_requires_sample_and_set(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--sample", "S1"],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
