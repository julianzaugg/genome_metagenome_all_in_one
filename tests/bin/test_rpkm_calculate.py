#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "bin" / "rpkm_calculate.py"
HEADER = "sample\tqseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tpercent_query_aligned\tpercent_subject_aligned\n"


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class RpkmCalculateTest(unittest.TestCase):
    def test_singlem_normalisation_and_zero_filled_gene_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            singlem_dir = tmp / "singlem"
            singlem_dir.mkdir()
            marker_lengths = tmp / "marker_lengths.tsv"
            gene_blast = tmp / "gene_blast.tsv"
            outdir = tmp / "out"

            marker_lengths.write_text(
                "Gene\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
                "S3.marker\t1\t300\t100\t100\t100\n"
            )
            (singlem_dir / "S3.marker_blast.tsv").write_text(
                HEADER
                + "S1\tread1\tmarker_hit\tmarker\t100\t50\t0\t0\t1\t50\t1\t50\t1e-20\t100\t150\t300\t0.333\t0.167\n"
                + "S1\tread2\tmarker_hit\tmarker\t100\t50\t0\t0\t1\t50\t1\t50\t1e-20\t100\t150\t300\t0.333\t0.167\n"
                + "S2\tread3\tmarker_hit\tmarker\t100\t50\t0\t0\t1\t50\t1\t50\t1e-20\t100\t150\t300\t0.333\t0.167\n"
            )
            gene_blast.write_text(
                HEADER
                + "S1\tread1\tgeneA\tgeneA\t100\t90\t0\t0\t1\t90\t1\t90\t1e-20\t100\t150\t300\t0.600\t0.300\n"
                + "S1\tread2\tgeneA\tgeneA\t100\t90\t0\t0\t1\t90\t1\t90\t1e-20\t100\t150\t300\t0.600\t0.300\n"
                + "S2\tread3\tgeneB\tgeneB\t100\t90\t0\t0\t1\t90\t1\t90\t1e-20\t100\t150\t300\t0.600\t0.300\n"
            )

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--singlem-blast-dir",
                    str(singlem_dir),
                    "--marker-lengths",
                    str(marker_lengths),
                    "--gene-blast",
                    str(gene_blast),
                    "--outdir",
                    str(outdir),
                ],
                check=True,
            )

            means = read_tsv(outdir / "singlem_rpkm_means.tsv")
            self.assertEqual({row["sample"] for row in means}, {"S1", "S2"})
            self.assertGreater(float(means[0]["Mean_rpkm"]), 0)

            counts = read_tsv(outdir / "gene_catalogue_mapped_reads_per_gene.tsv")
            by_gene = {row["Gene_ID"]: row for row in counts}
            self.assertEqual(by_gene["geneA"]["S1"], "2")
            self.assertEqual(by_gene["geneA"]["S2"], "0")
            self.assertEqual(by_gene["geneB"]["S1"], "0")
            self.assertEqual(by_gene["geneB"]["S2"], "1")

            normalised = read_tsv(outdir / "gene_catalogue_rpkm_per_gene_normalised.tsv")
            self.assertEqual({row["Gene_ID"] for row in normalised}, {"geneA", "geneB"})


if __name__ == "__main__":
    unittest.main()
