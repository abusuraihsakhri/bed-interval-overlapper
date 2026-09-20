"use strict";

const $ = (id) => document.getElementById(id);

const state = {
  worker: null,
  workerReady: false,
  initPromise: null,
  pending: new Map(),
  requestId: 0,
  lastResult: null,
  lastOperation: null,
};

const sampleA = `chrom,start,end,name,score,strand
chr1,10000,10500,PROM_BRCA1_01,100,+
chr1,10400,11000,PROM_BRCA1_02,95,+
chr1,25000,26000,ENH_MYC_01,80,-
chr2,15000,16500,EXON_TP53_01,100,+`;

const sampleB = `chrom,start,end,name,score,strand
chr1,10200,10800,H3K4ME3_01,90,+
chr1,25800,27000,PEAK_MYC_01,85,-
chr2,16400,17200,PEAK_TP53_01,90,+
chr3,100000,101500,CTCF_01,99,+`;

function setTheme(theme) {
  document.documentElement.dataset.theme = theme;
  localStorage.setItem("bed-theme", theme);
  $("theme-toggle").setAttribute("aria-label", theme === "dark" ? "Switch to light theme" : "Switch to dark theme");
}

function initializeTheme() {
  setTheme(localStorage.getItem("bed-theme") === "dark" ? "dark" : "light");
}

function setRuntimeStatus(text, ready = false) {
  const el = $("runtime-status");
  el.textContent = text;
  el.classList.toggle("ready", ready);
}

function setBusy(busy) {
  const button = $("analyze");
  button.disabled = busy;
  button.textContent = busy ? "Analyzing…" : "Analyze intervals";
}

function showMessage(message) {
  const el = $("message");
  el.textContent = message;
  el.hidden = !message;
}

function getWorker() {
  if (state.worker) return state.worker;
  state.worker = new Worker("worker.js");
  state.worker.addEventListener("message", (event) => {
    const { type, id, result, error } = event.data || {};
    if (type === "ready") {
      state.workerReady = true;
      setRuntimeStatus("Python runtime ready", true);
      return;
    }
    if (type === "status") {
      setRuntimeStatus(event.data.message || "Loading Python runtime…");
      return;
    }
    const pending = state.pending.get(id);
    if (!pending) return;
    state.pending.delete(id);
    if (error) pending.reject(new Error(error));
    else pending.resolve(result);
  });
  state.worker.addEventListener("error", (event) => {
    setRuntimeStatus("Runtime failed to load");
    showMessage(event.message || "The browser Python runtime could not be initialized.");
  });
  return state.worker;
}

function callWorker(payload) {
  const worker = getWorker();
  const id = ++state.requestId;
  return new Promise((resolve, reject) => {
    state.pending.set(id, { resolve, reject });
    worker.postMessage({ type: "analyze", id, payload });
  });
}

function requiresSetB(operation) {
  return ["intersect", "subtract", "jaccard"].includes(operation);
}

function updateOperationUI() {
  const operation = $("operation").value;
  const needsB = requiresSetB(operation);
  $("set-b-card").classList.toggle("hidden", !needsB);
  $("dataset-grid").classList.toggle("single", !needsB);
  $("intersect-options").classList.toggle("hidden", operation !== "intersect");
  $("merge-options").classList.toggle("hidden", operation !== "merge");
  $("coverage-options").classList.toggle("hidden", operation !== "coverage");
  updateInputCount();
}

function roughCount(text, format) {
  const lines = text.split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
  if (!lines.length) return 0;
  if (format === "csv" || (format === "auto" && lines[0].includes(","))) return Math.max(0, lines.length - 1);
  return lines.filter((line) => !line.startsWith("#") && !line.startsWith("track") && !line.startsWith("browser")).length;
}

function updateInputCount() {
  const format = $("input-format").value;
  const a = roughCount($("set-a").value, format);
  const b = requiresSetB($("operation").value) ? roughCount($("set-b").value, format) : 0;
  $("input-note").textContent = requiresSetB($("operation").value) ? `Approx. ${a} + ${b} intervals` : `Approx. ${a} intervals`;
}

function numericValue(id, fallback) {
  const value = Number($(id).value);
  return Number.isFinite(value) ? value : fallback;
}

function buildPayload() {
  return {
    operation: $("operation").value,
    format: $("input-format").value,
    setA: $("set-a").value,
    setB: $("set-b").value,
    options: {
      minOverlap: numericValue("min-overlap", 1),
      fractionA: numericValue("fraction-a", 0),
      fractionB: numericValue("fraction-b", 0),
      reciprocal: $("reciprocal").checked,
      strandMode: $("strand-mode").value,
      mergeDistance: numericValue("merge-distance", 0),
      binSize: numericValue("bin-size", 1000),
    },
  };
}

function validatePayload(payload) {
  if (!payload.setA.trim()) throw new Error("Set A is empty.");
  if (requiresSetB(payload.operation) && !payload.setB.trim()) throw new Error("Set B is required for this operation.");
  if (payload.options.minOverlap < 1 || !Number.isInteger(payload.options.minOverlap)) throw new Error("Minimum overlap must be an integer of at least 1 bp.");
  for (const [name, value] of [["Fraction A", payload.options.fractionA], ["Fraction B", payload.options.fractionB]]) {
    if (value < 0 || value > 1) throw new Error(`${name} must be between 0 and 1.`);
  }
  if (payload.options.binSize < 1 || !Number.isInteger(payload.options.binSize)) throw new Error("Bin size must be a positive integer.");
}

function formatNumber(value) {
  if (typeof value !== "number") return String(value ?? "—");
  return Number.isInteger(value) ? value.toLocaleString() : value.toLocaleString(undefined, { maximumFractionDigits: 6 });
}

function setSummary(items) {
  const cards = $("summary-cards");
  cards.innerHTML = "";
  for (const [label, value] of items.slice(0, 4)) {
    const card = document.createElement("div");
    card.className = "metric-card";
    const span = document.createElement("span");
    span.textContent = label;
    const strong = document.createElement("strong");
    strong.textContent = formatNumber(value);
    card.append(span, strong);
    cards.appendChild(card);
  }
}

function renderTable(columns, rows) {
  const table = $("result-table");
  const thead = table.querySelector("thead");
  const tbody = table.querySelector("tbody");
  thead.innerHTML = "";
  tbody.innerHTML = "";
  $("empty-state").hidden = rows.length > 0;
  table.hidden = rows.length === 0;
  if (!rows.length) return;

  const headerRow = document.createElement("tr");
  columns.forEach(([key, label]) => {
    const th = document.createElement("th");
    th.scope = "col";
    th.textContent = label;
    headerRow.appendChild(th);
  });
  thead.appendChild(headerRow);

  for (const row of rows.slice(0, 500)) {
    const tr = document.createElement("tr");
    columns.forEach(([key]) => {
      const td = document.createElement("td");
      td.textContent = formatNumber(row[key]);
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  }
}

function renderResult(result, operation) {
  state.lastResult = result;
  state.lastOperation = operation;
  $("download-results").disabled = false;
  showMessage("");

  if (operation === "intersect") {
    const rows = result.rows || [];
    setSummary([["Status", "Complete"], ["Set A", result.set_a_count], ["Set B", result.set_b_count], ["Overlaps", rows.length]]);
    renderTable([
      ["chrom", "Chrom"], ["overlap_start", "Start"], ["overlap_end", "End"], ["overlap_bp", "Overlap bp"],
      ["name_a", "A name"], ["name_b", "B name"], ["fraction_a", "Frac A"], ["fraction_b", "Frac B"],
    ], rows);
  } else if (operation === "merge") {
    const rows = result.rows || [];
    setSummary([["Status", "Complete"], ["Input", result.input_count], ["Merged", rows.length], ["Covered bp", result.total_bp]]);
    renderTable([["chrom", "Chrom"], ["start", "Start"], ["end", "End"], ["length", "Length"], ["name", "Name"]], rows);
  } else if (operation === "subtract") {
    const rows = result.rows || [];
    setSummary([["Status", "Complete"], ["Set A", result.set_a_count], ["Set B", result.set_b_count], ["Segments", rows.length]]);
    renderTable([["chrom", "Chrom"], ["start", "Start"], ["end", "End"], ["length", "Length"], ["name", "Name"], ["strand", "Strand"]], rows);
  } else if (operation === "jaccard") {
    const j = result.jaccard || {};
    setSummary([["Status", "Complete"], ["Jaccard", j.jaccard_index], ["Intersection bp", j.intersection_bp], ["Union bp", j.union_bp]]);
    renderTable([["metric", "Metric"], ["value", "Value"]], [
      { metric: "Set A intervals", value: j.set_a_intervals },
      { metric: "Set B intervals", value: j.set_b_intervals },
      { metric: "Set A merged bp", value: j.set_a_merged_bp },
      { metric: "Set B merged bp", value: j.set_b_merged_bp },
      { metric: "Overlap count", value: j.overlap_count },
    ]);
  } else if (operation === "coverage") {
    const rows = result.rows || [];
    setSummary([["Status", "Complete"], ["Intervals", result.total_intervals], ["Covered bp", result.total_covered_bp], ["Chromosomes", rows.length]]);
    renderTable([
      ["chrom", "Chrom"], ["span_start", "Span start"], ["span_end", "Span end"], ["span_bp", "Span bp"],
      ["covered_bp", "Covered bp"], ["breadth_coverage_fraction", "Breadth"], ["mean_depth", "Mean depth"], ["max_depth", "Max depth"],
    ], rows);
  }
}

async function analyze() {
  setBusy(true);
  showMessage("");
  try {
    const payload = buildPayload();
    validatePayload(payload);
    if (!state.workerReady) setRuntimeStatus("Loading Python runtime…");
    const result = await callWorker(payload);
    renderResult(result, payload.operation);
  } catch (error) {
    showMessage(error?.message || String(error));
    setSummary([["Status", "Error"], ["Set A", 0], ["Set B", 0], ["Result", "—"]]);
  } finally {
    setBusy(false);
  }
}

function clearAll() {
  $("set-a").value = "";
  $("set-b").value = "";
  state.lastResult = null;
  state.lastOperation = null;
  $("download-results").disabled = true;
  $("result-table").hidden = true;
  $("result-table").querySelector("thead").innerHTML = "";
  $("result-table").querySelector("tbody").innerHTML = "";
  $("empty-state").hidden = false;
  setSummary([["Status", "Ready"], ["Set A", 0], ["Set B", 0], ["Result", "—"]]);
  showMessage("");
  updateInputCount();
}

function loadSample() {
  $("set-a").value = sampleA;
  $("set-b").value = sampleB;
  $("input-format").value = "auto";
  updateInputCount();
}

async function loadFile(which) {
  const input = which === "a" ? $("file-a") : $("file-b");
  const file = input.files?.[0];
  if (!file) return;
  const text = await file.text();
  $(which === "a" ? "set-a" : "set-b").value = text;
  updateInputCount();
  input.value = "";
}

function downloadResults() {
  if (!state.lastResult) return;
  const payload = JSON.stringify(state.lastResult, null, 2);
  const blob = new Blob([payload], { type: "application/json;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `bed-${state.lastOperation || "result"}.json`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function bindEvents() {
  $("theme-toggle").addEventListener("click", () => setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark"));
  $("operation").addEventListener("change", updateOperationUI);
  $("input-format").addEventListener("change", updateInputCount);
  $("set-a").addEventListener("input", updateInputCount);
  $("set-b").addEventListener("input", updateInputCount);
  $("analyze").addEventListener("click", analyze);
  $("clear").addEventListener("click", clearAll);
  $("load-sample").addEventListener("click", loadSample);
  $("download-results").addEventListener("click", downloadResults);
  document.querySelectorAll("[data-file-target]").forEach((button) => {
    button.addEventListener("click", () => $(button.dataset.fileTarget === "a" ? "file-a" : "file-b").click());
  });
  $("file-a").addEventListener("change", () => loadFile("a"));
  $("file-b").addEventListener("change", () => loadFile("b"));
}

initializeTheme();
bindEvents();
updateOperationUI();
$("result-table").hidden = true;
