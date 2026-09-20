"use strict";

const PYODIDE_VERSION = "314.0.7";
const PYODIDE_BASE = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
let pyodidePromise = null;

function postStatus(message) {
  self.postMessage({ type: "status", message });
}

async function initialize() {
  if (pyodidePromise) return pyodidePromise;
  pyodidePromise = (async () => {
    postStatus("Loading Python runtime…");
    importScripts(`${PYODIDE_BASE}pyodide.js`);
    const pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });
    postStatus("Loading interval engine…");
    const response = await fetch("bed_overlap.py", { cache: "no-cache" });
    if (!response.ok) throw new Error(`Could not load bed_overlap.py (${response.status}).`);
    pyodide.FS.writeFile("bed_overlap.py", await response.text());
    await pyodide.runPythonAsync("import sys; sys.path.insert(0, '.'); import bed_overlap");
    self.postMessage({ type: "ready" });
    return pyodide;
  })();
  return pyodidePromise;
}

async function runAnalysis(payload) {
  const pyodide = await initialize();
  pyodide.globals.set("payload_json", JSON.stringify(payload));
  const output = await pyodide.runPythonAsync(`
import json
from bed_overlap import BedParser, IntervalEngine

payload = json.loads(payload_json)
fmt = payload.get("format", "auto")

def parse_set(text, slot):
    if fmt == "bed":
        return BedParser.parse_bed_string(text)
    path = f"/tmp/bed_interval_{slot}.csv"
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return BedParser.load_intervals(path)

set_a = parse_set(payload.get("setA", ""), "a")
set_b = parse_set(payload.get("setB", ""), "b") if payload.get("setB", "").strip() else []
if not set_a:
    raise ValueError("Set A contains no valid intervals.")

operation = payload["operation"]
options = payload.get("options", {})

if operation == "intersect":
    if not set_b:
        raise ValueError("Set B contains no valid intervals.")
    overlaps = IntervalEngine.intersect(
        set_a,
        set_b,
        min_overlap_bp=int(options.get("minOverlap", 1)),
        fraction_a=float(options.get("fractionA", 0.0)),
        fraction_b=float(options.get("fractionB", 0.0)),
        reciprocal=bool(options.get("reciprocal", False)),
        strand_mode=options.get("strandMode", "any"),
    )
    result = {
        "set_a_count": len(set_a),
        "set_b_count": len(set_b),
        "rows": [{
            "chrom": item["chrom"],
            "overlap_start": item["overlap_start"],
            "overlap_end": item["overlap_end"],
            "overlap_bp": item["overlap_bp"],
            "name_a": item["interval_a"]["name"],
            "name_b": item["interval_b"]["name"],
            "fraction_a": item["fraction_a"],
            "fraction_b": item["fraction_b"],
        } for item in overlaps],
    }
elif operation == "merge":
    merged = IntervalEngine.merge(set_a, max_distance=int(options.get("mergeDistance", 0)))
    rows = [item.to_dict() for item in merged]
    result = {"input_count": len(set_a), "total_bp": sum(item["length"] for item in rows), "rows": rows}
elif operation == "subtract":
    if not set_b:
        raise ValueError("Set B contains no valid intervals.")
    rows = [item.to_dict() for item in IntervalEngine.subtract(set_a, set_b)]
    result = {"set_a_count": len(set_a), "set_b_count": len(set_b), "rows": rows}
elif operation == "jaccard":
    if not set_b:
        raise ValueError("Set B contains no valid intervals.")
    result = {"jaccard": IntervalEngine.jaccard_similarity(set_a, set_b)}
elif operation == "coverage":
    coverage = IntervalEngine.coverage_profile(set_a, bin_size=int(options.get("binSize", 1000)))
    rows = [{"chrom": chrom, **values} for chrom, values in coverage["chromosomes"].items()]
    result = {"total_intervals": coverage["total_intervals"], "total_covered_bp": coverage["total_covered_bp"], "rows": rows}
else:
    raise ValueError(f"Unsupported operation: {operation}")

json.dumps(result)
`);
  return JSON.parse(output);
}

self.addEventListener("message", async (event) => {
  const { type, id, payload } = event.data || {};
  if (type !== "analyze") return;
  try {
    const result = await runAnalysis(payload);
    self.postMessage({ id, result });
  } catch (error) {
    self.postMessage({ id, error: error?.message || String(error) });
  }
});
