# BED Interval Overlapper & Genomic Operations Toolkit

A high-performance, pure Python genomic interval calculation engine implementing sweep-line interval operations, reciprocal overlap filters, coverage profiling, and Jaccard similarity metrics conforming to the **UCSC BED Format Specification** (0-based half-open coordinates $[start, end)$).

---

## Genomic Coordinate Conventions & Mathematical Models

### 1. Zero-Based Half-Open Interval Algebra
A genomic feature spanning chromosome coordinates from start to end represents nucleotides:
$$\text{Coordinates: } [start, end) \implies \text{Length} = end - start$$

Given two overlapping intervals $A = [s_A, e_A)$ and $B = [s_B, e_B)$ on the same chromosome:
$$\text{Overlap Length: } \text{len}(A \cap B) = \max(0, \min(e_A, e_B) - \max(s_A, s_B))$$
$$\text{Overlap Fraction } f_A = \frac{\text{len}(A \cap B)}{\text{len}(A)}, \quad f_B = \frac{\text{len}(A \cap B)}{\text{len}(B)}$$

### 2. Reciprocal Overlap & Jaccard Similarity
- **Reciprocal Overlap Constraint:** Validated if $f_A \ge \theta$ and $f_B \ge \theta$ where $\theta \in (0, 1]$.
- **Genomic Jaccard Similarity Index:**
$$J(A_{\text{set}}, B_{\text{set}}) = \frac{\text{Total Intersecting Base Pairs}}{\text{Total Union Base Pairs}} = \frac{\sum |A_i \cap B_j|}{\sum |A_i \cup B_j|}$$

### 3. Sweep-Line Interval Merging & Coverage
Sorted interval boundaries are processed as event points $E = \{(s_i, +1), (e_i, -1)\}$. Active interval count increments at starts and decrements at ends, computing depth profile bins and merging intervals with distance gap tolerance $d$.

---

## Features

- **Standard Operations:** Full suite equivalent to `bedtools`:
  - `intersect`: Intersection with minimum overlap, fraction $f$, and reciprocal constraints.
  - `merge`: Linear-time sweep-line merging with user-defined gap distance.
  - `subtract`: Complementary subtraction yielding remaining fragments.
  - `jaccard`: Base-pair intersection-over-union metric.
  - `coverage`: Binned depth profiling across genomic tracks.
- **Dual Format Ingestion:** Seamlessly reads standard `.bed` format (BED3 to BED6) and tabular `.csv` files.
- **High-Throughput Batch Processing:** Command-line batch execution for genomic feature annotation pipelines.
- **Zero Runtime Dependencies:** Standalone implementation utilizing only the Python Standard Library.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies.

```bash
git clone https://github.com/abusuraihsakhri/bed-interval-overlapper.git
cd bed-interval-overlapper
```

---

## CLI Usage

### 1. Merge Overlapping Genomic Intervals
```bash
python cli.py merge -i sample.csv --distance 50
```

### 2. Intersect Two Interval Sets
```bash
python cli.py intersect -a sample.csv -b target.bed --min-overlap 100 --reciprocal 0.5
```

### 3. Calculate Genomic Jaccard Index
```bash
python cli.py jaccard -a set1.bed -b set2.bed
```

### 4. Batch Process CSV/BED Files
```bash
python cli.py batch --input sample.csv --output results.csv
```

---

## Python API Quickstart

```python
from bed_overlap import IntervalEngine, BedParser

intervals = BedParser.parse_csv_file("sample.csv")
merged = IntervalEngine.merge(intervals, max_distance=0)

print(f"Original: {len(intervals)} intervals -> Merged: {len(merged)} intervals")
for m in merged:
    print(f"  {m.chrom}:{m.start}-{m.end} (len={len(m)})")
```

---

## Testing & Verification

Run the test suite:

```bash
python -m pytest -p no:zarr
```

