#!/usr/bin/env python3
"""
BED Interval Overlapper & Genomic Interval Operations Toolkit.

Implements:
- 0-based half-open [start, end) BED interval parsing (BED3 to BED6+)
- Sweep-line interval intersection (equivalent to bedtools intersect)
- Reciprocal overlap fraction, min overlap length, and strand-aware matching
- Interval merging with gap tolerance (equivalent to bedtools merge)
- Interval subtraction / difference (equivalent to bedtools subtract)
- Jaccard Similarity Index calculation for genomic feature sets
- Binned coverage depth profiling & breadth of coverage statistics
- BEDPE paired-end fragment spanning and concordance filtering
- SVG genome browser track visualizer

Pure Python Standard Library (no external dependencies required).
"""

from __future__ import annotations
import math
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple, Optional, Any, Set


@dataclass(frozen=True)
class BedInterval:
    """Genomic interval in 0-based half-open coordinate space [start, end)."""
    chrom: str
    start: int
    end: int
    name: str = "."
    score: str = "."
    strand: str = "."  # '+', '-', or '.'

    def __post_init__(self):
        if self.start < 0:
            raise ValueError(f"Interval start ({self.start}) cannot be negative.")
        if self.end < self.start:
            raise ValueError(f"Interval end ({self.end}) cannot be less than start ({self.start}).")

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlaps_with(self, other: BedInterval) -> bool:
        """Check if two intervals overlap on the same chromosome."""
        if self.chrom != other.chrom:
            return False
        return max(self.start, other.start) < min(self.end, other.end)

    def overlap_bp(self, other: BedInterval) -> int:
        """Calculate overlap size in base pairs."""
        if self.chrom != other.chrom:
            return 0
        return max(0, min(self.end, other.end) - max(self.start, other.start))

    def to_bed_line(self) -> str:
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.name}\t{self.score}\t{self.strand}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "name": self.name,
            "score": self.score,
            "strand": self.strand,
        }


@dataclass(frozen=True)
class BedpeRecord:
    """Paired-end genomic interval (BEDPE)."""
    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: str = "."
    score: str = "."
    strand1: str = "+"
    strand2: str = "-"

    @property
    def is_intra_chromosomal(self) -> bool:
        return self.chrom1 == self.chrom2

    @property
    def insert_size(self) -> int:
        if not self.is_intra_chromosomal:
            return -1
        return max(self.end1, self.end2) - min(self.start1, self.start2)

    def spanning_interval(self) -> BedInterval:
        if not self.is_intra_chromosomal:
            chrom_label = f"{self.chrom1}|{self.chrom2}"
        else:
            chrom_label = self.chrom1
        s = min(self.start1, self.start2)
        e = max(self.end1, self.end2)
        return BedInterval(chrom_label, s, e, name=self.name, score=self.score)


class BedParser:
    """Parses and formats BED format files and strings."""

    @staticmethod
    def parse_line(line: str) -> Optional[BedInterval]:
        """Parse a single BED line into BedInterval."""
        clean = line.strip()
        if not clean or clean.startswith("#") or clean.startswith("track") or clean.startswith("browser"):
            return None
        fields = clean.split("\t")
        if len(fields) < 3:
            fields = clean.split()
        if len(fields) < 3:
            return None

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3] if len(fields) > 3 else "."
        score = fields[4] if len(fields) > 4 else "."
        strand = fields[5] if len(fields) > 5 else "."
        if strand not in ("+", "-", "."):
            strand = "."

        return BedInterval(chrom, start, end, name, score, strand)

    @classmethod
    def parse_bed_string(cls, content: str) -> List[BedInterval]:
        """Parse multi-line BED text into list of BedIntervals."""
        intervals = []
        for line in content.strip().splitlines():
            iv = cls.parse_line(line)
            if iv:
                intervals.append(iv)
        return intervals

    @classmethod
    def load_bed_file(cls, filepath: str) -> List[BedInterval]:
        """Read BED file from disk."""
        intervals = []
        with open(filepath, "r", encoding="utf-8-sig") as f:
            for line in f:
                iv = cls.parse_line(line)
                if iv:
                    intervals.append(iv)
        return intervals


class IntervalEngine:
    """High-performance genomic interval operations."""

    @staticmethod
    def group_by_chromosome(intervals: List[BedInterval]) -> Dict[str, List[BedInterval]]:
        """Group and sort intervals by chromosome."""
        grouped: Dict[str, List[BedInterval]] = {}
        for iv in intervals:
            grouped.setdefault(iv.chrom, []).append(iv)
        for chrom in grouped:
            grouped[chrom].sort(key=lambda x: (x.start, x.end))
        return grouped

    @classmethod
    def intersect(
        cls,
        set_a: List[BedInterval],
        set_b: List[BedInterval],
        min_overlap_bp: int = 1,
        fraction_a: float = 0.0,
        fraction_b: float = 0.0,
        reciprocal: bool = False,
        strand_mode: str = "any",  # 'any', 'same', 'opposite'
    ) -> List[Dict[str, Any]]:
        """
        Sweep-line intersection between two interval sets.
        Supports min overlap bp, fractional overlap of A/B, reciprocal fraction, and strandedness.
        """
        grouped_a = cls.group_by_chromosome(set_a)
        grouped_b = cls.group_by_chromosome(set_b)
        results: List[Dict[str, Any]] = []

        for chrom, list_a in grouped_a.items():
            if chrom not in grouped_b:
                continue
            list_b = grouped_b[chrom]

            j = 0
            len_b = len(list_b)

            for a in list_a:
                len_a = a.length
                if len_a == 0:
                    continue

                # Advance pointer j until list_b[j].end > a.start
                while j < len_b and list_b[j].end <= a.start:
                    j += 1

                k = j
                while k < len_b and list_b[k].start < a.end:
                    b = list_b[k]
                    k += 1

                    # Strand check
                    if strand_mode == "same" and a.strand != "." and b.strand != "." and a.strand != b.strand:
                        continue
                    elif strand_mode == "opposite" and a.strand != "." and b.strand != "." and a.strand == b.strand:
                        continue

                    # Overlap span
                    ov_start = max(a.start, b.start)
                    ov_end = min(a.end, b.end)
                    ov_bp = max(0, ov_end - ov_start)

                    if ov_bp < min_overlap_bp:
                        continue

                    frac_a = ov_bp / len_a
                    len_b_val = b.length or 1
                    frac_b = ov_bp / len_b_val

                    if fraction_a > 0.0 and frac_a < fraction_a:
                        continue
                    if fraction_b > 0.0 and frac_b < fraction_b:
                        continue
                    if reciprocal:
                        req = max(fraction_a, fraction_b)
                        if frac_a < req or frac_b < req:
                            continue

                    results.append({
                        "chrom": chrom,
                        "overlap_start": ov_start,
                        "overlap_end": ov_end,
                        "overlap_bp": ov_bp,
                        "fraction_a": round(frac_a, 4),
                        "fraction_b": round(frac_b, 4),
                        "interval_a": a.to_dict(),
                        "interval_b": b.to_dict(),
                    })

        return results

    @classmethod
    def merge(cls, intervals: List[BedInterval], max_distance: int = 0) -> List[BedInterval]:
        """
        Merge overlapping or adjacent intervals within max_distance bp.
        Equivalent to bedtools merge.
        """
        grouped = cls.group_by_chromosome(intervals)
        merged_all: List[BedInterval] = []

        for chrom in sorted(grouped.keys()):
            iv_list = grouped[chrom]
            if not iv_list:
                continue

            current_start = iv_list[0].start
            current_end = iv_list[0].end
            merged_names = [iv_list[0].name] if iv_list[0].name != "." else []

            for iv in iv_list[1:]:
                if iv.start <= (current_end + max_distance):
                    current_end = max(current_end, iv.end)
                    if iv.name != "." and iv.name not in merged_names:
                        merged_names.append(iv.name)
                else:
                    name_str = ";".join(merged_names) if merged_names else "."
                    merged_all.append(BedInterval(chrom, current_start, current_end, name=name_str))
                    current_start = iv.start
                    current_end = iv.end
                    merged_names = [iv.name] if iv.name != "." else []

            name_str = ";".join(merged_names) if merged_names else "."
            merged_all.append(BedInterval(chrom, current_start, current_end, name=name_str))

        return merged_all

    @classmethod
    def subtract(cls, set_a: List[BedInterval], set_b: List[BedInterval]) -> List[BedInterval]:
        """
        Subtract intervals in set_b from set_a.
        Equivalent to bedtools subtract.
        """
        merged_b = cls.merge(set_b)
        grouped_b = cls.group_by_chromosome(merged_b)
        result: List[BedInterval] = []

        for a in set_a:
            if a.chrom not in grouped_b:
                result.append(a)
                continue

            b_list = grouped_b[a.chrom]
            current_segments = [(a.start, a.end)]

            for b in b_list:
                if b.end <= a.start or b.start >= a.end:
                    continue

                next_segments = []
                for s, e in current_segments:
                    # No overlap with current piece
                    if b.end <= s or b.start >= e:
                        next_segments.append((s, e))
                    else:
                        # Left remainder
                        if s < b.start:
                            next_segments.append((s, b.start))
                        # Right remainder
                        if e > b.end:
                            next_segments.append((b.end, e))
                current_segments = next_segments
                if not current_segments:
                    break

            for s, e in current_segments:
                if s < e:
                    result.append(BedInterval(a.chrom, s, e, name=a.name, score=a.score, strand=a.strand))

        return result

    @classmethod
    def jaccard_similarity(cls, set_a: List[BedInterval], set_b: List[BedInterval]) -> Dict[str, Any]:
        """
        Calculate Jaccard similarity index:
        J(A, B) = Intersection_bp / Union_bp
        where Union_bp = Total_bp(Merge(A)) + Total_bp(Merge(B)) - Intersection_bp
        """
        merged_a = cls.merge(set_a)
        merged_b = cls.merge(set_b)

        bp_a = sum(iv.length for iv in merged_a)
        bp_b = sum(iv.length for iv in merged_b)

        overlaps = cls.intersect(merged_a, merged_b)
        intersection_bp = sum(o["overlap_bp"] for o in overlaps)

        union_bp = bp_a + bp_b - intersection_bp
        jaccard = round(intersection_bp / union_bp, 6) if union_bp > 0 else 0.0

        return {
            "set_a_intervals": len(set_a),
            "set_b_intervals": len(set_b),
            "set_a_merged_bp": bp_a,
            "set_b_merged_bp": bp_b,
            "intersection_bp": intersection_bp,
            "union_bp": union_bp,
            "jaccard_index": jaccard,
            "overlap_count": len(overlaps),
        }

    @classmethod
    def coverage_profile(cls, intervals: List[BedInterval], bin_size: int = 1000) -> Dict[str, Any]:
        """
        Calculate genomic coverage depth profile using sweep-line event points.
        """
        if not intervals:
            return {"total_intervals": 0, "total_covered_bp": 0, "chromosomes": {}}

        grouped = cls.group_by_chromosome(intervals)
        chrom_profiles: Dict[str, Any] = {}
        total_covered_bp = 0

        for chrom, iv_list in grouped.items():
            events: List[Tuple[int, int]] = []
            for iv in iv_list:
                events.append((iv.start, 1))
                events.append((iv.end, -1))
            events.sort(key=lambda x: (x[0], -x[1]))

            min_coord = min(iv.start for iv in iv_list)
            max_coord = max(iv.end for iv in iv_list)
            span_bp = max_coord - min_coord

            depth = 0
            cur_pos = events[0][0]
            covered_bases = 0
            max_depth = 0
            depth_times_length = 0

            binned_depths: Dict[int, int] = {}

            for pos, delta in events:
                if pos > cur_pos and depth > 0:
                    length = pos - cur_pos
                    covered_bases += length
                    depth_times_length += depth * length
                    max_depth = max(max_depth, depth)

                    b0 = cur_pos // bin_size
                    b1 = (pos - 1) // bin_size
                    for b in range(b0, b1 + 1):
                        bs = b * bin_size
                        be = (b + 1) * bin_size
                        w = min(be, pos) - max(bs, cur_pos)
                        if w > 0:
                            binned_depths[b] = binned_depths.get(b, 0) + depth * w

                depth += delta
                cur_pos = pos

            total_covered_bp += covered_bases
            mean_depth = round(depth_times_length / span_bp, 3) if span_bp > 0 else 0.0
            breadth_cov = round(covered_bases / span_bp, 4) if span_bp > 0 else 0.0

            chrom_profiles[chrom] = {
                "span_start": min_coord,
                "span_end": max_coord,
                "span_bp": span_bp,
                "covered_bp": covered_bases,
                "breadth_coverage_fraction": breadth_cov,
                "mean_depth": mean_depth,
                "max_depth": max_depth,
                "bins_count": len(binned_depths),
            }

        return {
            "total_intervals": len(intervals),
            "total_covered_bp": total_covered_bp,
            "chromosomes": chrom_profiles,
        }

    @staticmethod
    def render_svg(
        set_a: List[BedInterval],
        set_b: List[BedInterval],
        output_svg_path: str,
        width: int = 1000,
        height: int = 400
    ) -> None:
        """Render SVG track visualization of interval sets A and B."""
        all_chroms = sorted(set(list(IntervalEngine.group_by_chromosome(set_a).keys()) +
                                list(IntervalEngine.group_by_chromosome(set_b).keys())))

        svg = ET.Element("svg", xmlns="http://www.w3.org/2000/svg", width=str(width), height=str(height))
        style = ET.SubElement(svg, "style")
        style.text = """
            .title { font: bold 14px sans-serif; fill: #1e293b; }
            .track-label { font: bold 12px sans-serif; fill: #334155; }
            .ruler { stroke: #cbd5e1; stroke-width: 1; }
            .ruler-text { font: 10px monospace; fill: #64748b; }
            .iv-a { fill: #3b82f6; opacity: 0.85; }
            .iv-b { fill: #ef4444; opacity: 0.85; }
        """

        ET.SubElement(svg, "text", x="20", y="30", class_="title").text = "Genomic BED Interval Overlap Tracks"

        y = 60
        margin = 120
        track_w = width - margin - 40

        for chrom in all_chroms[:4]:
            a_list = [iv for iv in set_a if iv.chrom == chrom]
            b_list = [iv for iv in set_b if iv.chrom == chrom]

            all_starts = [iv.start for iv in a_list + b_list]
            all_ends = [iv.end for iv in a_list + b_list]
            if not all_starts:
                continue

            c_min = min(all_starts)
            c_max = max(all_ends)
            span = max(1, c_max - c_min)
            scale = track_w / span

            # Chromosome label & ruler
            ET.SubElement(svg, "text", x="20", y=str(y + 14), class_="track-label").text = f"{chrom}"
            ET.SubElement(svg, "line", x1=str(margin), y1=str(y + 8), x2=str(margin + track_w), y2=str(y + 8), class_="ruler")

            # Track A
            y_a = y + 18
            for iv in a_list:
                x = margin + int((iv.start - c_min) * scale)
                w = max(2, int(iv.length * scale))
                ET.SubElement(svg, "rect", x=str(x), y=str(y_a), width=str(w), height="10", rx="2", class_="iv-a")

            # Track B
            y_b = y + 34
            for iv in b_list:
                x = margin + int((iv.start - c_min) * scale)
                w = max(2, int(iv.length * scale))
                ET.SubElement(svg, "rect", x=str(x), y=str(y_b), width=str(w), height="10", rx="2", class_="iv-b")

            y += 65

        tree = ET.ElementTree(svg)
        ET.indent(tree, space="  ")
        tree.write(output_svg_path, encoding="unicode", xml_declaration=True)


def get_curated_bed_benchmarks() -> Dict[str, Dict[str, List[BedInterval]]]:
    """Curated benchmark BED sets (ChIP-seq peaks, promoter regions, enhancer elements)."""
    return {
        "promoters_vs_h3k4me3": {
            "set_a": [  # Promoter regions
                BedInterval("chr1", 10000, 11500, "PROM_01", "100", "+"),
                BedInterval("chr1", 25000, 26200, "PROM_02", "100", "-"),
                BedInterval("chr1", 50000, 51000, "PROM_03", "100", "+"),
                BedInterval("chr2", 15000, 16500, "PROM_04", "100", "+"),
                BedInterval("chr2", 80000, 81500, "PROM_05", "100", "-"),
            ],
            "set_b": [  # H3K4me3 ChIP-seq active promoter peaks
                BedInterval("chr1", 10500, 12000, "PEAK_01", "90", "+"),
                BedInterval("chr1", 24800, 26000, "PEAK_02", "95", "-"),
                BedInterval("chr1", 90000, 91500, "PEAK_03", "70", "+"),
                BedInterval("chr2", 15200, 16800, "PEAK_04", "85", "+"),
            ]
        },
        "ctcf_insulators_vs_super_enhancers": {
            "set_a": [  # CTCF binding sites
                BedInterval("chr3", 100000, 100500, "CTCF_01", "100", "+"),
                BedInterval("chr3", 200000, 200500, "CTCF_02", "100", "-"),
                BedInterval("chr3", 300000, 300500, "CTCF_03", "100", "+"),
            ],
            "set_b": [  # Super enhancers
                BedInterval("chr3", 100200, 105000, "SE_01", "90", "."),
                BedInterval("chr3", 299000, 302000, "SE_02", "80", "."),
            ]
        }
    }
