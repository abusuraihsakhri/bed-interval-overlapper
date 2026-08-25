#!/usr/bin/env python3
"""
Unit Test Suite for BED Interval Overlapper & Genomic Arithmetic Engine.
Covers:
- Coordinate validity & interval constraints [start, end)
- BedParser line and file format parsing (BED3 to BED6)
- Sweep-line interval intersection with min overlap, fractional, reciprocal, and stranded filters
- Interval merging (adjacent, overlapping, gap distance parameters)
- Interval subtraction and multiple segment splitting
- Jaccard similarity index calculation (identical, disjoint, partial)
- Genomic coverage depth profiling (breadth, mean depth, max depth)
- BEDPE paired-end record spanning
- SVG track rendering
- Curated genomic benchmarks
- CLI subcommands, batch workflows, and JSON serialization
"""

import json
import os
import sys
import unittest
from io import StringIO
from unittest.mock import patch

from bed_overlap import (
    BedInterval,
    BedpeRecord,
    BedParser,
    IntervalEngine,
    get_curated_bed_benchmarks,
)
import cli


class TestBedIntervalAndParser(unittest.TestCase):
    def test_valid_interval(self):
        iv = BedInterval("chr1", 100, 200, "GENE1", "100", "+")
        self.assertEqual(iv.chrom, "chr1")
        self.assertEqual(iv.start, 100)
        self.assertEqual(iv.end, 200)
        self.assertEqual(iv.length, 100)
        self.assertEqual(iv.strand, "+")

    def test_invalid_interval_bounds(self):
        with self.assertRaises(ValueError):
            BedInterval("chr1", -10, 100)
        with self.assertRaises(ValueError):
            BedInterval("chr1", 200, 100)

    def test_parse_bed_lines(self):
        # BED3
        iv3 = BedParser.parse_line("chr1\t1000\t2000")
        self.assertIsNotNone(iv3)
        self.assertEqual(iv3.chrom, "chr1")
        self.assertEqual(iv3.start, 1000)
        self.assertEqual(iv3.end, 2000)

        # BED6
        iv6 = BedParser.parse_line("chr2\t5000\t6000\tPEAK1\t99\t-")
        self.assertIsNotNone(iv6)
        self.assertEqual(iv6.name, "PEAK1")
        self.assertEqual(iv6.strand, "-")

        # Comment and header lines should return None
        self.assertIsNone(BedParser.parse_line("# Header line"))
        self.assertIsNone(BedParser.parse_line("track name=test"))
        self.assertIsNone(BedParser.parse_line("browser position chr1:1-1000"))
        self.assertIsNone(BedParser.parse_line(""))


class TestIntervalIntersection(unittest.TestCase):
    def setUp(self):
        self.set_a = [
            BedInterval("chr1", 100, 300, "A1", ".", "+"),
            BedInterval("chr1", 500, 700, "A2", ".", "-"),
            BedInterval("chr2", 1000, 2000, "A3", ".", "+"),
        ]
        self.set_b = [
            BedInterval("chr1", 200, 400, "B1", ".", "+"),
            BedInterval("chr1", 550, 650, "B2", ".", "+"),
            BedInterval("chr2", 3000, 4000, "B3", ".", "+"),
        ]

    def test_basic_intersection(self):
        overlaps = IntervalEngine.intersect(self.set_a, self.set_b)
        self.assertEqual(len(overlaps), 2)

        # A1 (100-300) with B1 (200-400) -> overlap 200-300 (100 bp)
        ov1 = next(o for o in overlaps if o["interval_a"]["name"] == "A1")
        self.assertEqual(ov1["overlap_start"], 200)
        self.assertEqual(ov1["overlap_end"], 300)
        self.assertEqual(ov1["overlap_bp"], 100)

        # A2 (500-700) with B2 (550-650) -> overlap 550-650 (100 bp)
        ov2 = next(o for o in overlaps if o["interval_a"]["name"] == "A2")
        self.assertEqual(ov2["overlap_start"], 550)
        self.assertEqual(ov2["overlap_end"], 650)
        self.assertEqual(ov2["overlap_bp"], 100)

    def test_min_overlap_filter(self):
        overlaps = IntervalEngine.intersect(self.set_a, self.set_b, min_overlap_bp=150)
        self.assertEqual(len(overlaps), 0)

    def test_strand_matching_same(self):
        # A2 is '-', B2 is '+' -> should be excluded under strand_mode="same"
        overlaps = IntervalEngine.intersect(self.set_a, self.set_b, strand_mode="same")
        self.assertEqual(len(overlaps), 1)
        self.assertEqual(overlaps[0]["interval_a"]["name"], "A1")

    def test_strand_matching_opposite(self):
        # A2 is '-', B2 is '+' -> included under strand_mode="opposite"
        overlaps = IntervalEngine.intersect(self.set_a, self.set_b, strand_mode="opposite")
        self.assertEqual(len(overlaps), 1)
        self.assertEqual(overlaps[0]["interval_a"]["name"], "A2")

    def test_fractional_overlap_and_reciprocal(self):
        # A1 len=200, B1 len=200, overlap=100 -> frac_a = 0.5, frac_b = 0.5
        ov_rec = IntervalEngine.intersect(self.set_a, self.set_b, fraction_a=0.5, fraction_b=0.5, reciprocal=True)
        self.assertTrue(any(o["interval_a"]["name"] == "A1" for o in ov_rec))

        ov_high = IntervalEngine.intersect(self.set_a, self.set_b, fraction_a=0.9, reciprocal=True)
        self.assertEqual(len(ov_high), 0)


class TestIntervalMergeAndSubtract(unittest.TestCase):
    def test_merge_overlapping_intervals(self):
        ivs = [
            BedInterval("chr1", 100, 200, "N1"),
            BedInterval("chr1", 150, 300, "N2"),
            BedInterval("chr1", 400, 500, "N3"),
            BedInterval("chr1", 500, 600, "N4"),  # Adjacent
            BedInterval("chr2", 50, 100, "N5"),
        ]
        merged = IntervalEngine.merge(ivs)
        self.assertEqual(len(merged), 3)

        # chr1: 100-300 (merged N1, N2)
        self.assertEqual(merged[0].chrom, "chr1")
        self.assertEqual(merged[0].start, 100)
        self.assertEqual(merged[0].end, 300)

        # chr1: 400-600 (merged N3, N4)
        self.assertEqual(merged[1].chrom, "chr1")
        self.assertEqual(merged[1].start, 400)
        self.assertEqual(merged[1].end, 600)

        # chr2: 50-100
        self.assertEqual(merged[2].chrom, "chr2")
        self.assertEqual(merged[2].start, 50)
        self.assertEqual(merged[2].end, 100)

    def test_merge_with_distance_gap(self):
        ivs = [
            BedInterval("chr1", 100, 200),
            BedInterval("chr1", 250, 350),  # Gap of 50 bp
        ]
        # max_distance = 60 should merge them into 100-350
        merged = IntervalEngine.merge(ivs, max_distance=60)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].start, 100)
        self.assertEqual(merged[0].end, 350)

    def test_subtract_middle_split(self):
        # Set A: [100, 500)
        # Set B: [200, 300)
        # Result should be [100, 200) and [300, 500)
        a = [BedInterval("chr1", 100, 500, "GENE")]
        b = [BedInterval("chr1", 200, 300, "EXON")]
        sub = IntervalEngine.subtract(a, b)
        self.assertEqual(len(sub), 2)
        self.assertEqual(sub[0].start, 100)
        self.assertEqual(sub[0].end, 200)
        self.assertEqual(sub[1].start, 300)
        self.assertEqual(sub[1].end, 500)

    def test_subtract_complete_overlap(self):
        a = [BedInterval("chr1", 100, 200)]
        b = [BedInterval("chr1", 50, 250)]
        sub = IntervalEngine.subtract(a, b)
        self.assertEqual(len(sub), 0)


class TestJaccardAndCoverage(unittest.TestCase):
    def test_jaccard_identical_sets(self):
        a = [BedInterval("chr1", 100, 200), BedInterval("chr1", 300, 400)]
        b = [BedInterval("chr1", 100, 200), BedInterval("chr1", 300, 400)]
        jac = IntervalEngine.jaccard_similarity(a, b)
        self.assertEqual(jac["jaccard_index"], 1.0)
        self.assertEqual(jac["intersection_bp"], 200)
        self.assertEqual(jac["union_bp"], 200)

    def test_jaccard_disjoint_sets(self):
        a = [BedInterval("chr1", 100, 200)]
        b = [BedInterval("chr1", 300, 400)]
        jac = IntervalEngine.jaccard_similarity(a, b)
        self.assertEqual(jac["jaccard_index"], 0.0)
        self.assertEqual(jac["intersection_bp"], 0)

    def test_coverage_depth_profile(self):
        ivs = [
            BedInterval("chr1", 100, 300),  # len 200
            BedInterval("chr1", 200, 400),  # len 200, overlap 200-300 (depth 2)
        ]
        cov = IntervalEngine.coverage_profile(ivs, bin_size=100)
        self.assertEqual(cov["total_intervals"], 2)
        self.assertEqual(cov["total_covered_bp"], 300)  # span 100 to 400 = 300 bp
        self.assertEqual(cov["chromosomes"]["chr1"]["max_depth"], 2)
        self.assertEqual(cov["chromosomes"]["chr1"]["breadth_coverage_fraction"], 1.0)


class TestBedpeAndVisualization(unittest.TestCase):
    def test_bedpe_record(self):
        pe = BedpeRecord("chr1", 100, 150, "chr1", 400, 450, "PAIR1")
        self.assertTrue(pe.is_intra_chromosomal)
        self.assertEqual(pe.insert_size, 350)
        span = pe.spanning_interval()
        self.assertEqual(span.chrom, "chr1")
        self.assertEqual(span.start, 100)
        self.assertEqual(span.end, 450)

    def test_svg_rendering(self):
        a = [BedInterval("chr1", 1000, 2000, "A")]
        b = [BedInterval("chr1", 1500, 2500, "B")]
        tmp_svg = "test_output_track.svg"
        try:
            IntervalEngine.render_svg(a, b, tmp_svg)
            self.assertTrue(os.path.exists(tmp_svg))
            with open(tmp_svg, "r") as f:
                content = f.read()
            self.assertIn("<svg", content)
            self.assertIn("chr1", content)
        finally:
            if os.path.exists(tmp_svg):
                os.remove(tmp_svg)


class TestCLIInterface(unittest.TestCase):
    def setUp(self):
        self.bed_a = "test_a.bed"
        self.bed_b = "test_b.bed"
        with open(self.bed_a, "w") as f:
            f.write("chr1\t1000\t2000\tA1\t100\t+\n")
            f.write("chr1\t3000\t4000\tA2\t100\t-\n")
        with open(self.bed_b, "w") as f:
            f.write("chr1\t1500\t2500\tB1\t100\t+\n")
            f.write("chr1\t5000\t6000\tB2\t100\t+\n")

    def tearDown(self):
        for p in (self.bed_a, self.bed_b):
            if os.path.exists(p):
                os.remove(p)

    def test_cli_intersect(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["intersect", "--a", self.bed_a, "--b", self.bed_b])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("chr1\t1500\t2000\tA1|B1\t500", output)

    def test_cli_intersect_json(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["intersect", "--a", self.bed_a, "--b", self.bed_b, "--json"])
            self.assertEqual(ret, 0)
            data = json.loads(fake_out.getvalue())
            self.assertIsInstance(data, list)
            self.assertEqual(len(data), 1)
            self.assertEqual(data[0]["overlap_bp"], 500)

    def test_cli_merge(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["merge", "--input", self.bed_a])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("chr1\t1000\t2000", output)

    def test_cli_subtract(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["subtract", "--a", self.bed_a, "--b", self.bed_b])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("chr1\t1000\t1500", output)

    def test_cli_jaccard(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["jaccard", "--a", self.bed_a, "--b", self.bed_b])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("Jaccard Index:", output)

    def test_cli_coverage(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["coverage", "--input", self.bed_a])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("Genomic Coverage Profile", output)

    def test_cli_benchmark(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["benchmark", "promoters_vs_h3k4me3"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("GENOMIC INTERVAL BENCHMARK REPORT", output)

    def test_cli_list_benchmarks(self):
        with patch('sys.stdout', new=StringIO()) as fake_out:
            ret = cli.main(["--list-benchmarks"])
            self.assertEqual(ret, 0)
            output = fake_out.getvalue()
            self.assertIn("promoters_vs_h3k4me3", output)


if __name__ == "__main__":
    unittest.main()
