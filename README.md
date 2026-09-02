# Bed Interval Overlapper

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

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

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`BedInterval`**: Genomic interval in 0-based half-open coordinate space [start, end).
- **`BedpeRecord`**: Paired-end genomic interval (BEDPE).
- **`BedParser`**: Parses and formats BED format files and strings.
- **`IntervalEngine`**: High-performance genomic interval operations.

---

## 📐 Mathematical Formulation & Logic

```text
  """Calculate overlap size in base pairs."""
  score = fields[4] if len(fields) > 4 else "."
  Calculate Jaccard similarity index:
  Calculate genomic coverage depth profile using sweep-line event points.
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --a <value> --b <value> --min-overlap <value> --input <value>
```

### Parameter Reference
- `--a`: Specifies input measurement or parameter value.
- `--b`: Specifies input measurement or parameter value.
- `--min-overlap`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--distance`: Specifies input measurement or parameter value.
- `--bin-size`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.
- `--interactive`: Specifies input measurement or parameter value.
- `---`: Specifies input measurement or parameter value.
- `--fraction-a`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `suite_name` | Parameter / observation metric | Required |
| `system_slug` | Parameter / observation metric | Required |
| `standard_reference` | Parameter / observation metric | Required |
| `test_cases` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t bed-interval-overlapper .
docker run -p 8000:8000 bed-interval-overlapper
```
