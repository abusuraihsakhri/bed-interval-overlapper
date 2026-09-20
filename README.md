# BED Interval Overlapper

A dependency-free Python toolkit for common genomic interval operations on 0-based, half-open BED coordinates. The repository includes a command-line interface and a browser interface that runs the same Python interval engine locally with Pyodide.

## Features

- Intersect two interval sets with minimum base-pair overlap, fractional overlap, reciprocal filtering, and strand matching.
- Merge overlapping or book-ended intervals with an optional distance tolerance.
- Subtract one interval set from another.
- Calculate base-pair Jaccard similarity.
- Calculate coverage breadth and depth summaries.
- Read BED3–BED6 and simple CSV interval tables.
- Render SVG interval tracks from the Python CLI.
- Run the core Python engine entirely in the browser; input data remains local to the browser.

## Browser application

Open `index.html` through a local web server or use the GitHub Pages deployment once enabled. The interface accepts pasted data or local BED/CSV files and supports light and dark themes.

The browser loads a pinned Pyodide runtime from jsDelivr. Interval contents are processed locally and are not uploaded by the application.

## CLI

Python 3.10+ is supported and the runtime has no third-party Python dependencies.

```bash
python cli.py merge -i sample.csv --distance 50
python cli.py intersect -a set_a.bed -b set_b.bed --min-overlap 100 --fraction-a 0.5 --fraction-b 0.5 --reciprocal
python cli.py jaccard -a set_a.bed -b set_b.bed
python cli.py coverage -i set_a.bed --bin-size 1000
python cli.py batch -i sample.csv -o results.csv
```

## Python API

```python
from bed_overlap import BedParser, IntervalEngine

intervals = BedParser.load_intervals("sample.csv")
merged = IntervalEngine.merge(intervals, max_distance=0)

for interval in merged:
    print(interval.chrom, interval.start, interval.end, interval.length)
```

BED coordinates follow the standard 0-based, half-open convention: `[start, end)` and `length = end - start`.

## Local development

```bash
python -m pip install pytest
python -m pytest -v
python -m http.server 8000
```

Then open `http://localhost:8000/` for the browser interface.

## Browser compatibility

The browser application requires a modern browser with WebAssembly, Web Workers, `fetch`, and ES2020-era JavaScript support. Pyodide is pinned to version 314.0.7.

## License

MIT. See [LICENSE](LICENSE).
