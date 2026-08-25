# BED Interval Overlapper & Genomic Arithmetic Engine

A zero-dependency, high-performance genomic arithmetic and interval manipulation engine implementing sweep-line interval intersections, gap-tolerant interval merging, interval subtractions, Jaccard similarity index calculation, depth of coverage profiling, and SVG genome browser track rendering.

## Core Operations

### 1. Interval Intersection (`bedtools intersect` equivalent)
- **Coordinate System**: Standard 0-based half-open coordinates $[start, end)$.
- **Sweep-Line Overlap Detection**: Fast $O(N \log N)$ chromosome-partitioned sweep-line traversal.
- **Filtering Capabilities**:
  - Minimum overlap length in base pairs (`--min-overlap`).
  - Fractional overlap threshold of query A (`--fraction-a`) or target B (`--fraction-b`).
  - Reciprocal overlap requirement (`--reciprocal`).
  - Strandedness matching: `any` (ignore strand), `same` (same strand only), `opposite` (opposite strand only).

### 2. Interval Merging (`bedtools merge` equivalent)
- Merges overlapping or book-ended ($end_1 = start_2$) intervals on the same chromosome.
- Supports maximum distance gap merging (`--distance <bp>`).

### 3. Interval Subtraction (`bedtools subtract` equivalent)
- Subtracts overlapping regions in set B from intervals in set A, splitting remaining non-overlapping segments.

### 4. Jaccard Similarity Index (`bedtools jaccard` equivalent)
$$J(A, B) = \frac{|\text{Intersection}(A, B)|}{|\text{Union}(A, B)|} = \frac{\text{Overlap BP}}{\text{Merged BP}(A) + \text{Merged BP}(B) - \text{Overlap BP}}$$

### 5. Coverage Depth & Breadth Profiling (`genomecov` equivalent)
- Computes base-pair resolution and binned depth counts, breadth of coverage fraction, mean sequencing depth, and peak maximum depth.

### 6. BEDPE (Paired-End BED) & SVG Track Visualizer
- Handles paired-end sequencing fragments and intra-chromosomal spans.
- Generates stand-alone SVG track diagrams representing comparative genomic features.

---

## Benchmark Datasets

Includes built-in curated genomic feature comparisons:
1. **promoters_vs_h3k4me3**: Promoter annotations vs active histone mark ChIP-seq peak calls.
2. **ctcf_insulators_vs_super_enhancers**: CTCF boundary insulators vs super-enhancer genomic spans.

---

## CLI Usage

### 1. Intersect BED Files
```bash
python cli.py intersect --a promoters.bed --b peaks.bed --min-overlap 50
```

### 2. Merge Intervals with Gap Distance
```bash
python cli.py merge --input peaks.bed --distance 100
```

### 3. Subtract Target Intervals
```bash
python cli.py subtract --a genes.bed --b repeats.bed
```

### 4. Jaccard Similarity
```bash
python cli.py jaccard --a sample_a.bed --b sample_b.bed
```

### 5. Coverage Depth Profiling
```bash
python cli.py coverage --input peaks.bed --bin-size 1000
```

### 6. SVG Genome Track Rendering
```bash
python cli.py visualize --a set_a.bed --b set_b.bed --output tracks.svg
```

### 7. Run Built-in Benchmark
```bash
python cli.py benchmark promoters_vs_h3k4me3
```

### 8. Interactive Shell
```bash
python cli.py --interactive
```

---

## Unit Testing

Run unit tests via Python's standard `unittest` framework:

```bash
python -m unittest test_bed_overlap.py
```

All 25 test cases validate boundary conditions, sweep-line correctness, reciprocal overlap calculations, multi-segment subtractions, Jaccard indices, and CLI commands.
