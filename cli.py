#!/usr/bin/env python3
"""
Command Line Interface for BED Interval Overlapper & Genomic Interval Engine.

Usage examples:
  python cli.py intersect --a promoters.bed --b peaks.bed --min-overlap 50
  python cli.py merge --input intervals.bed --distance 100
  python cli.py subtract --a genes.bed --b repeats.bed
  python cli.py jaccard --a sample_a.bed --b sample_b.bed
  python cli.py coverage --input peaks.bed --bin-size 500
  python cli.py benchmark promoters_vs_h3k4me3
  python cli.py benchmark promoters_vs_h3k4me3 --json
  python cli.py --interactive
"""

import argparse
import json
import os
import sys
from typing import List, Optional

from bed_overlap import (
    BedInterval,
    BedParser,
    IntervalEngine,
    get_curated_bed_benchmarks,
)


def format_intervals_table(intervals: List[BedInterval]) -> str:
    lines = ["  Chrom       Start         End       Length  Name                 Strand"]
    lines.append("  " + "-" * 75)
    for iv in intervals[:50]:
        lines.append(f"  {iv.chrom:10s} {iv.start:>10,d} {iv.end:>11,d} {iv.length:>10,d}  {iv.name:20s} {iv.strand}")
    if len(intervals) > 50:
        lines.append(f"  ... and {len(intervals) - 50} more intervals.")
    return "\n".join(lines)


def run_interactive():
    benchmarks = get_curated_bed_benchmarks()
    print("\n" + "=" * 80)
    print("  BED INTERVAL OVERLAPPER & GENOMIC ARITHMETIC - INTERACTIVE SHELL")
    print("=" * 80)
    print("Type 'help' for commands, 'exit' or 'quit' to exit.\n")

    while True:
        try:
            cmd_line = input("bed-tools> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExiting.")
            break

        if not cmd_line:
            continue

        parts = cmd_line.split()
        cmd = parts[0].lower()

        if cmd in ("exit", "quit", "q"):
            print("Goodbye.")
            break
        elif cmd == "help":
            print("""
Available Commands:
  benchmarks                      List curated benchmark dataset pairs
  eval-bm <benchmark_id>          Run intersect, jaccard & coverage on a benchmark pair
  help                            Show this menu
  exit / quit                     Exit the shell
            """)
        elif cmd == "benchmarks":
            print("\nAvailable Curated Genomic Benchmark Pairs:")
            for k, v in benchmarks.items():
                print(f"  [{k:35s}] Set A: {len(v['set_a'])} intervals | Set B: {len(v['set_b'])} intervals")
            print()
        elif cmd == "eval-bm":
            if len(parts) < 2:
                print("Usage: eval-bm <benchmark_id>")
                continue
            b_id = parts[1]
            if b_id not in benchmarks:
                print(f"Benchmark '{b_id}' not found. Use 'benchmarks' to see list.")
                continue
            bm = benchmarks[b_id]
            set_a, set_b = bm["set_a"], bm["set_b"]
            overlaps = IntervalEngine.intersect(set_a, set_b)
            jaccard = IntervalEngine.jaccard_similarity(set_a, set_b)
            cov_a = IntervalEngine.coverage_profile(set_a)
            print(f"\n--- Benchmark Evaluation: {b_id} ---")
            print(f"Set A Count: {len(set_a)} | Set B Count: {len(set_b)}")
            print(f"Overlapping Regions: {len(overlaps)}")
            print(f"Jaccard Similarity:  {jaccard['jaccard_index']:.6f} (Intersection: {jaccard['intersection_bp']:,} bp / Union: {jaccard['union_bp']:,} bp)")
            print(f"Set A Covered Bases: {cov_a['total_covered_bp']:,} bp\n")
        else:
            print(f"Unknown command '{cmd}'. Type 'help' for instructions.")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="bed-interval-overlapper",
        description="BED Interval Overlapper & Genomic Arithmetic Engine",
    )

    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # Intersect
    p_inter = subparsers.add_parser("intersect", help="Find overlaps between BED set A and set B")
    p_inter.add_argument("--a", required=True, help="Path to BED file A (query)")
    p_inter.add_argument("--b", required=True, help="Path to BED file B (target)")
    p_inter.add_argument("--min-overlap", type=int, default=1, help="Minimum overlap in bp (default: 1)")
    p_inter.add_argument("--fraction-a", type=float, default=0.0, help="Minimum fraction of A that must overlap")
    p_inter.add_argument("--fraction-b", type=float, default=0.0, help="Minimum fraction of B that must overlap")
    p_inter.add_argument("--reciprocal", action="store_true", help="Require reciprocal overlap fraction")
    p_inter.add_argument("--strand", choices=["any", "same", "opposite"], default="any", help="Strand matching mode")
    p_inter.add_argument("--json", action="store_true", help="Output in JSON format")

    # Merge
    p_merge = subparsers.add_parser("merge", help="Merge overlapping/adjacent intervals in a BED file")
    p_merge.add_argument("--input", "-i", required=True, help="Path to BED file to merge")
    p_merge.add_argument("--distance", "-d", type=int, default=0, help="Maximum distance gap to merge across (default: 0)")
    p_merge.add_argument("--json", action="store_true", help="Output in JSON format")

    # Subtract
    p_sub = subparsers.add_parser("subtract", help="Subtract intervals in set B from set A")
    p_sub.add_argument("--a", required=True, help="Path to BED file A")
    p_sub.add_argument("--b", required=True, help="Path to BED file B to subtract")
    p_sub.add_argument("--json", action="store_true", help="Output in JSON format")

    # Jaccard
    p_jac = subparsers.add_parser("jaccard", help="Calculate Jaccard similarity index between set A and B")
    p_jac.add_argument("--a", required=True, help="Path to BED file A")
    p_jac.add_argument("--b", required=True, help="Path to BED file B")
    p_jac.add_argument("--json", action="store_true", help="Output in JSON format")

    # Coverage
    p_cov = subparsers.add_parser("coverage", help="Calculate genomic coverage depth profile")
    p_cov.add_argument("--input", "-i", required=True, help="Path to BED file")
    p_cov.add_argument("--bin-size", type=int, default=1000, help="Bin size for depth profiling in bp")
    p_cov.add_argument("--json", action="store_true", help="Output in JSON format")

    # Visualize
    p_viz = subparsers.add_parser("visualize", help="Render SVG track visualization")
    p_viz.add_argument("--a", required=True, help="Path to BED file A")
    p_viz.add_argument("--b", required=True, help="Path to BED file B")
    p_viz.add_argument("--output", "-o", default="genome_tracks.svg", help="Output SVG path")

    # Benchmark
    p_bm = subparsers.add_parser("benchmark", help="Run operations on curated genomic benchmark datasets")
    p_bm.add_argument("benchmark_name", nargs="?", default="promoters_vs_h3k4me3", help="Curated benchmark pair name")
    p_bm.add_argument("--json", action="store_true", help="Output in JSON format")

    # Global flags
    parser.add_argument("--list-benchmarks", action="store_true", help="List all available curated benchmark pairs")
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive shell")

    args = parser.parse_args(argv)
    benchmarks = get_curated_bed_benchmarks()

    if args.interactive:
        run_interactive()
        return 0

    if args.list_benchmarks:
        if getattr(args, "json", False):
            print(json.dumps({k: {"set_a_count": len(v["set_a"]), "set_b_count": len(v["set_b"])}
                             for k, v in benchmarks.items()}, indent=2))
        else:
            print("\nAvailable Curated Genomic Benchmark Pairs:")
            print("=" * 75)
            for k, v in benchmarks.items():
                print(f"  [{k:35s}] Set A: {len(v['set_a'])} intervals | Set B: {len(v['set_b'])} intervals")
            print("=" * 75)
        return 0

    if args.subcommand == "benchmark":
        bm_name = args.benchmark_name
        if bm_name not in benchmarks:
            err = {"error": f"Benchmark '{bm_name}' not found. Choose from: {list(benchmarks.keys())}"}
            if args.json:
                print(json.dumps(err, indent=2))
            else:
                print(f"Error: {err['error']}", file=sys.stderr)
            return 1

        bm = benchmarks[bm_name]
        set_a, set_b = bm["set_a"], bm["set_b"]
        overlaps = IntervalEngine.intersect(set_a, set_b)
        jaccard = IntervalEngine.jaccard_similarity(set_a, set_b)
        cov_a = IntervalEngine.coverage_profile(set_a)
        cov_b = IntervalEngine.coverage_profile(set_b)

        res = {
            "benchmark": bm_name,
            "set_a_count": len(set_a),
            "set_b_count": len(set_b),
            "overlap_count": len(overlaps),
            "jaccard_similarity": jaccard,
            "set_a_coverage": cov_a,
            "set_b_coverage": cov_b,
            "overlaps": overlaps,
        }

        if args.json:
            print(json.dumps(res, indent=2))
        else:
            print("=" * 80)
            print(f"  GENOMIC INTERVAL BENCHMARK REPORT: [{bm_name}]")
            print("=" * 80)
            print(f"  Set A Intervals: {len(set_a)} | Set B Intervals: {len(set_b)}")
            print(f"  Total Overlaps:  {len(overlaps)}")
            print(f"  Jaccard Index:   {jaccard['jaccard_index']:.6f}")
            print(f"  Intersection:    {jaccard['intersection_bp']:,} bp")
            print(f"  Union:           {jaccard['union_bp']:,} bp")
            print("-" * 80)
            print("  Overlapping Intervals (First 5):")
            for ov in overlaps[:5]:
                print(f"    * {ov['chrom']}:{ov['overlap_start']}-{ov['overlap_end']} ({ov['overlap_bp']} bp) "
                      f"[A: {ov['interval_a']['name']} <==> B: {ov['interval_b']['name']}]")
            print("=" * 80)
        return 0

    if args.subcommand == "intersect":
        set_a = BedParser.load_bed_file(args.a)
        set_b = BedParser.load_bed_file(args.b)
        overlaps = IntervalEngine.intersect(
            set_a, set_b,
            min_overlap_bp=args.min_overlap,
            fraction_a=args.fraction_a,
            fraction_b=args.fraction_b,
            reciprocal=args.reciprocal,
            strand_mode=args.strand,
        )
        if args.json:
            print(json.dumps(overlaps, indent=2))
        else:
            print(f"\nIntersection Results ({len(overlaps)} overlaps found):")
            for ov in overlaps:
                print(f"{ov['chrom']}\t{ov['overlap_start']}\t{ov['overlap_end']}\t"
                      f"{ov['interval_a']['name']}|{ov['interval_b']['name']}\t{ov['overlap_bp']}\t.")
        return 0

    if args.subcommand == "merge":
        intervals = BedParser.load_bed_file(args.input)
        merged = IntervalEngine.merge(intervals, max_distance=args.distance)
        if args.json:
            print(json.dumps([iv.to_dict() for iv in merged], indent=2))
        else:
            print(f"\nMerged BED Intervals ({len(intervals)} -> {len(merged)}):")
            for iv in merged:
                print(iv.to_bed_line())
        return 0

    if args.subcommand == "subtract":
        set_a = BedParser.load_bed_file(args.a)
        set_b = BedParser.load_bed_file(args.b)
        subtracted = IntervalEngine.subtract(set_a, set_b)
        if args.json:
            print(json.dumps([iv.to_dict() for iv in subtracted], indent=2))
        else:
            print(f"\nSubtracted BED Intervals ({len(subtracted)} remaining segments):")
            for iv in subtracted:
                print(iv.to_bed_line())
        return 0

    if args.subcommand == "jaccard":
        set_a = BedParser.load_bed_file(args.a)
        set_b = BedParser.load_bed_file(args.b)
        jac = IntervalEngine.jaccard_similarity(set_a, set_b)
        if args.json:
            print(json.dumps(jac, indent=2))
        else:
            print(f"\nJaccard Similarity Index:")
            print(f"  Jaccard Index:   {jac['jaccard_index']:.6f}")
            print(f"  Intersection:    {jac['intersection_bp']:,} bp")
            print(f"  Union:           {jac['union_bp']:,} bp")
            print(f"  Overlap Count:   {jac['overlap_count']}")
        return 0

    if args.subcommand == "coverage":
        intervals = BedParser.load_bed_file(args.input)
        cov = IntervalEngine.coverage_profile(intervals, bin_size=args.bin_size)
        if args.json:
            print(json.dumps(cov, indent=2))
        else:
            print(f"\nGenomic Coverage Profile (Total Covered: {cov['total_covered_bp']:,} bp):")
            for chrom, data in cov["chromosomes"].items():
                print(f"  {chrom:10s} : Span: {data['span_bp']:,} bp | Covered: {data['covered_bp']:,} bp ({data['breadth_coverage_fraction']*100:.1f}%) | Mean Depth: {data['mean_depth']}x | Max Depth: {data['max_depth']}x")
        return 0

    if args.subcommand == "visualize":
        set_a = BedParser.load_bed_file(args.a)
        set_b = BedParser.load_bed_file(args.b)
        IntervalEngine.render_svg(set_a, set_b, args.output)
        print(f"SVG track visualization written to {args.output}")
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
