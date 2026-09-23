#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "build_read_stat_report.py"


def read_tsv(path):
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return rows[0], rows[1:]


def write(path, text):
    path.write_text(text)


class MetagenomeReportTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "seqkit").mkdir()
        (self.root / "repmags").mkdir()
        (self.root / "hq").mkdir()
        (self.root / "assess").mkdir()
        (self.root / "assess_genomes").mkdir()

        write(self.root / "seqkit" / "JZ_01.raw.seqkit_stats.tsv",
              "file\tformat\ttype\tnum_seqs\tsum_len\n"
              "JZ_01_R1.fastq.gz\tFASTQ\tDNA\t1000\t150000\n"
              "JZ_01_R2.fastq.gz\tFASTQ\tDNA\t1000\t150000\n")

        write(self.root / "repmags" / "JZ_01_abundances.tsv",
              "Genome\tJZ_01 Relative Abundance (%)\tJZ_01 Covered Fraction\tJZ_01 Mean\tJZ_01 Count\n"
              "unmapped\t10.0\tNA\tNA\tNA\n"
              "mag1\t60.0\t0.9\t5.0\t500\n"
              "mag2\t30.0\t0.5\t2.0\t300\n")

        write(self.root / "hq" / "mag1.fasta", ">contig1\nACGT\n")

        write(self.root / "assess" / "JZ_01.Dereplicated_Bins.mapping_assessment.tsv",
              "sample\tset\treads_mapped\treads_total\tpct_reads_mapped\treads_supplementary\t"
              "bases_mapped_readlen\tpct_bases_mapped_A\tbases_mapped_cigar\tpct_bases_mapped_B\t"
              "bases_mapped_cigar_all\tpct_bases_mapped_C\tbases_total\ttotal_source\tnote\n"
              "JZ_01\tDereplicated_Bins\t800\t2000\t40.00\t5\t120000\t80.00\t110000\t73.33\t"
              "115000\t76.67\t150000\tsupplied\t\n")

        write(self.root / "assess_genomes" / "JZ_01.Dereplicated_Bins.mapping_assessment_per_genome.tsv",
              "set\tsample\tgenome\treads_primary\treads_supplementary\tbases_A\tpct_bases_mapped_A\t"
              "bases_B\tpct_bases_mapped_B\tbases_C\tpct_bases_mapped_C\n"
              "Dereplicated_Bins\tJZ_01\tmag1\t500\t3\t72000\t48.00\t66000\t44.00\t70000\t46.67\n"
              "Dereplicated_Bins\tJZ_01\tmag2\t300\t2\t48000\t32.00\t44000\t29.33\t45000\t30.00\n")

    def run_report(self, out_name="report.tsv", extra_args=()):
        out = self.root / out_name
        cmd = [sys.executable, str(SCRIPT), "--mode", "metagenome", "--out", str(out),
               "--seqkit-dir", str(self.root / "seqkit"),
               "--repmag-dir", str(self.root / "repmags"),
               "--hq-dir", str(self.root / "hq"),
               "--assess-dir", str(self.root / "assess"),
               "--assess-genome-dir", str(self.root / "assess_genomes"),
               *extra_args]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return read_tsv(out)

    def test_bases_mapped_columns_follow_reads_mapped_pair(self):
        header, rows = self.run_report()
        i = header.index("Reads_mapped_Dereplicated_Bins_count")
        self.assertEqual(header[i:i + 8], [
            "Reads_mapped_Dereplicated_Bins_count",
            "Reads_mapped_Dereplicated_Bins_percent",
            "Bases_mapped_A_Dereplicated_Bins_count",
            "Bases_mapped_A_Dereplicated_Bins_percent",
            "Bases_mapped_B_Dereplicated_Bins_count",
            "Bases_mapped_B_Dereplicated_Bins_percent",
            "Bases_mapped_C_Dereplicated_Bins_count",
            "Bases_mapped_C_Dereplicated_Bins_percent",
        ])

    def test_bases_mapped_values_use_raw_bases_denominator(self):
        header, rows = self.run_report()
        row = dict(zip(header, rows[0]))
        # raw bases = 300000 (2x150000); GBbp should reflect that, not raw reads.
        self.assertEqual(row["GBbp"], "0.00030")
        self.assertEqual(row["Bases_mapped_A_Dereplicated_Bins_count"], "120000")
        self.assertEqual(row["Bases_mapped_A_Dereplicated_Bins_percent"], "40.00")
        self.assertEqual(row["Bases_mapped_B_Dereplicated_Bins_count"], "110000")
        self.assertEqual(row["Bases_mapped_B_Dereplicated_Bins_percent"], "36.67")
        self.assertEqual(row["Bases_mapped_C_Dereplicated_Bins_count"], "115000")
        self.assertEqual(row["Bases_mapped_C_Dereplicated_Bins_percent"], "38.33")

    def test_hq_subset_summed_from_per_genome_table(self):
        header, rows = self.run_report()
        row = dict(zip(header, rows[0]))
        # HQ_MAGs subset = mag1 only (the only genome in --hq-dir).
        self.assertEqual(row["Reads_mapped_HQ_MAGs_count"], "500")
        self.assertEqual(row["Bases_mapped_A_HQ_MAGs_count"], "72000")
        self.assertEqual(row["Bases_mapped_A_HQ_MAGs_percent"], "24.00")
        self.assertEqual(row["Bases_mapped_B_HQ_MAGs_count"], "66000")
        self.assertEqual(row["Bases_mapped_C_HQ_MAGs_count"], "70000")

    def test_bases_metrics_flag_trims_columns(self):
        header, rows = self.run_report("report_ab.tsv", extra_args=["--bases-metrics", "A,B"])
        bases_cols = [h for h in header if h.startswith("Bases_mapped_")]
        self.assertTrue(all("_C_" not in h for h in bases_cols))
        self.assertTrue(any("_A_" in h for h in bases_cols))
        self.assertTrue(any("_B_" in h for h in bases_cols))

    def test_invalid_bases_metric_errors(self):
        out = self.root / "bad.tsv"
        cmd = [sys.executable, str(SCRIPT), "--mode", "metagenome", "--out", str(out),
               "--seqkit-dir", str(self.root / "seqkit"), "--bases-metrics", "A,Z"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)

    def test_missing_assess_data_yields_blank_cells_not_crash(self):
        out = self.root / "no_assess.tsv"
        cmd = [sys.executable, str(SCRIPT), "--mode", "metagenome", "--out", str(out),
               "--seqkit-dir", str(self.root / "seqkit"),
               "--repmag-dir", str(self.root / "repmags"),
               "--hq-dir", str(self.root / "hq")]
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        header, rows = read_tsv(out)
        row = dict(zip(header, rows[0]))
        self.assertEqual(row["Bases_mapped_A_Dereplicated_Bins_count"], "")
        self.assertEqual(row["Bases_mapped_A_Dereplicated_Bins_percent"], "")
        # the read-count columns must still be populated
        self.assertEqual(row["Reads_mapped_Dereplicated_Bins_count"], "800")


class HostRemovalStageTest(unittest.TestCase):
    """Cleanifier (host removal) counts appear as their own QC stage after the others."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.seqkit = Path(self.tmp.name) / "seqkit"
        self.seqkit.mkdir()

    def stage(self, sample, stage, rows):
        text = "file\tformat\ttype\tnum_seqs\tsum_len\n"
        text += "".join(f"{sample}_{i}.fastq.gz\tFASTQ\tDNA\t{n}\t{n * 150}\n" for i, n in enumerate(rows, 1))
        write(self.seqkit / f"{sample}.{stage}.seqkit_stats.tsv", text)

    def run_report(self):
        out = Path(self.tmp.name) / "report.tsv"
        result = subprocess.run([sys.executable, str(SCRIPT), "--mode", "metagenome", "--out", str(out),
                                 "--seqkit-dir", str(self.seqkit)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        header, rows = read_tsv(out)
        return header, {row[0]: dict(zip(header, row)) for row in rows}

    def test_paired_cleanifier_counts_are_summed_across_mates(self):
        self.stage("JZ_01", "raw", [1000, 1000])
        self.stage("JZ_01", "fastp", [900, 900])
        self.stage("JZ_01", "cleanifier", [100, 100])
        header, rows = self.run_report()
        self.assertEqual(header, ["Sample_ID", "GBbp", "Raw_count", "Fastp_count", "Fastp_percent",
                                  "Cleanifier_count", "Cleanifier_percent"])
        self.assertEqual(rows["JZ_01"]["Cleanifier_count"], "200")
        self.assertEqual(rows["JZ_01"]["Cleanifier_percent"], "10.00")

    def test_long_read_cleanifier_follows_porechop_and_fastplong(self):
        self.stage("barcode01", "raw_long", [3037244])
        self.stage("barcode01", "porechop", [3027635])
        self.stage("barcode01", "fastplong", [2573447])
        self.stage("barcode01", "cleanifier", [26314])
        header, rows = self.run_report()
        self.assertEqual(header[3:], ["Porechop_count", "Porechop_percent", "Fastplong_count", "Fastplong_percent",
                                      "Cleanifier_count", "Cleanifier_percent"])
        self.assertEqual(rows["barcode01"]["Cleanifier_count"], "26314")
        self.assertEqual(rows["barcode01"]["Cleanifier_percent"], "0.87")

    def test_no_cleanifier_columns_without_host_removal(self):
        self.stage("barcode01", "raw_long", [1000])
        self.stage("barcode01", "fastplong", [800])
        header, _ = self.run_report()
        self.assertFalse(any(name.startswith("Cleanifier") for name in header))


class IsolateReportTest(unittest.TestCase):
    def test_isolate_mode_adds_bases_mapped_next_to_read_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "seqkit").mkdir()
            (root / "scaffolds").mkdir()
            (root / "assess").mkdir()

            write(root / "seqkit" / "ISO_1.raw.seqkit_stats.tsv",
                  "file\tformat\ttype\tnum_seqs\tsum_len\n"
                  "ISO_1_R1.fastq.gz\tFASTQ\tDNA\t500\t75000\n")

            write(root / "scaffolds" / "ISO_1_counts.tsv",
                  "Genome\tISO_1 Covered Fraction\tISO_1 Mean\tISO_1 Count\n"
                  "assembly\t0.95\t20.0\t450\n")

            write(root / "assess" / "ISO_1.Scaffolds.mapping_assessment.tsv",
                  "sample\tset\treads_mapped\treads_total\tpct_reads_mapped\treads_supplementary\t"
                  "bases_mapped_readlen\tpct_bases_mapped_A\tbases_mapped_cigar\tpct_bases_mapped_B\t"
                  "bases_mapped_cigar_all\tpct_bases_mapped_C\tbases_total\ttotal_source\tnote\n"
                  "ISO_1\tScaffolds\t450\t500\t90.00\t0\t60000\t80.00\t58000\t77.33\t"
                  "58000\t77.33\t75000\tsupplied\t\n")

            out = root / "report.tsv"
            cmd = [sys.executable, str(SCRIPT), "--mode", "isolate", "--out", str(out),
                   "--seqkit-dir", str(root / "seqkit"),
                   "--scaffold-dir", str(root / "scaffolds"),
                   "--assess-dir", str(root / "assess")]
            result = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)

            header, rows = read_tsv(out)
            row = dict(zip(header, rows[0]))
            i = header.index("Read_count")
            self.assertEqual(header[i:i + 8], [
                "Read_count", "Read_count_percent",
                "Bases_mapped_A_Scaffolds_count", "Bases_mapped_A_Scaffolds_percent",
                "Bases_mapped_B_Scaffolds_count", "Bases_mapped_B_Scaffolds_percent",
                "Bases_mapped_C_Scaffolds_count", "Bases_mapped_C_Scaffolds_percent",
            ])
            self.assertEqual(row["Bases_mapped_A_Scaffolds_count"], "60000")
            self.assertEqual(row["Bases_mapped_A_Scaffolds_percent"], "80.00")


if __name__ == "__main__":
    unittest.main()
