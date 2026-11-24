# The page for transport_view.jl: the leg's main vessels and tissue in front view, coloured by their iodine, and the iodine in
# the leg over time, the two in step. `write_view(path, D)` takes the data as a NamedTuple (transport_view.jl builds it); it
# needs nothing else, so made-up data can check the layout. One self-contained file, no external scripts or fonts.

using Base64, Printf

_num(x) = (y = Float64(x); !isfinite(y) ? "null" : y == round(y) && abs(y) < 1.0e15 ? string(Int(round(y))) : @sprintf("%.6g", y))
_arr(v; digits = 3) = "[" * join((_num(round(Float64(x); digits)) for x in v), ",") * "]"
_str(s) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "</" => "<\\/") * "\""

"The page's data as JSON."
function view_json(D)
    io = IOBuffer()
    seg(s) = string(
        "{\"x1\":", _arr(s.x1; digits = 1), ",\"z1\":", _arr(s.z1; digits = 1), ",\"x2\":", _arr(s.x2; digits = 1),
        ",\"z2\":", _arr(s.z2; digits = 1), ",\"r\":", _arr(s.r; digits = 3),
    )
    print(io, "{\"t\":", _arr(D.t; digits = 2))
    print(io, ",\"aif\":", _arr(D.aif; digits = 4), ",\"aif_dt\":", _num(D.aif_dt), ",\"aif_max\":", _num(D.aif_max))
    print(io, ",\"art\":", seg(D.art), ",\"tau\":", _arr(D.art.tau; digits = 2), "}")
    print(io, ",\"vein\":", seg(D.vein), ",\"q\":", _str(base64encode(D.vein.q)), ",\"max\":", _num(D.vein.max), "}")
    print(io, ",\"tis\":{\"x\":", _arr(D.tis.x; digits = 1), ",\"z\":", _arr(D.tis.z; digits = 1), ",\"ml\":", _arr(D.tis.ml; digits = 2))
    print(io, ",\"q\":", _str(base64encode(D.tis.q)), ",\"max\":", _num(D.tis.max), "}")
    print(io, ",\"sil\":{\"h\":", _num(D.sil.h), ",\"rows\":[", join((_arr(r; digits = 1) for r in D.sil.rows), ","), "]}")
    print(io, ",\"bounds\":{\"x0\":", _num(D.bounds.x0), ",\"x1\":", _num(D.bounds.x1), ",\"z0\":", _num(D.bounds.z0), ",\"z1\":", _num(D.bounds.z1), "}")
    print(io, ",\"series\":{\"keys\":[", join(_str.(D.series.keys), ","), "],\"names\":[", join(_str.(D.series.names), ","), "],\"values\":[")
    print(io, join((_arr(v; digits = 2) for v in D.series.values), ","), "]}")
    print(io, ",\"meta\":{", join((_str(String(k)) * ":" * (v isa AbstractString ? _str(v) : _num(v)) for (k, v) in pairs(D.meta)), ","), "}")
    print(io, "}")
    return String(take!(io))
end

const VIEW_TEMPLATE = raw"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Leg contrast at rest</title>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --leg: #efeee9; --empty: #d3d2ca;
  --vein: #2a78d6; --artery: #eb6834; --tissue: #1baf7a; --total: #52514e;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --leg: #242422; --empty: #3d3d3a;
    --vein: #3987e5; --artery: #d95926; --tissue: #199e70; --total: #c3c2b7;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --leg: #242422; --empty: #3d3d3a;
  --vein: #3987e5; --artery: #d95926; --tissue: #199e70; --total: #c3c2b7;
}
* { box-sizing: border-box; }
html, body { margin: 0; background: var(--page); color: var(--ink); }
body { font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1240px; margin: 0 auto; padding: 20px 16px 32px; }
h1 { font-size: 20px; font-weight: 600; margin: 0 0 4px; }
h2 { font-size: 14px; font-weight: 600; margin: 0 0 8px; }
.sub { color: var(--ink2); margin: 0 0 14px; max-width: 72ch; }
.controls { display: flex; flex-wrap: wrap; align-items: center; gap: 10px 18px; margin: 0 0 14px; }
.controls button { font: inherit; font-weight: 600; min-width: 84px; padding: 6px 14px; border-radius: 8px; border: 1px solid var(--border); background: var(--surface); color: var(--ink); cursor: pointer; }
.controls button:focus-visible, .controls input:focus-visible, .controls select:focus-visible, svg:focus-visible { outline: 2px solid var(--vein); outline-offset: 2px; }
.time { display: flex; align-items: center; gap: 10px; flex: 1 1 320px; min-width: 220px; }
.time input[type=range] { flex: 1; accent-color: var(--ink2); }
.readout { font-variant-numeric: tabular-nums; min-width: 8ch; color: var(--ink); font-weight: 600; }
.controls label { color: var(--ink2); display: inline-flex; align-items: center; gap: 6px; }
.controls select { font: inherit; padding: 4px 6px; border-radius: 6px; border: 1px solid var(--border); background: var(--surface); color: var(--ink); }
.panels { display: flex; flex-wrap: wrap; gap: 16px; align-items: flex-start; }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 14px 16px; }
.legcard { flex: 0 1 420px; min-width: 240px; }
.chartcard { flex: 1 1 520px; min-width: 280px; position: relative; }
.legwrap { display: flex; flex-wrap: wrap; gap: 12px 18px; align-items: flex-start; }
canvas { display: block; max-width: 100%; }
.scales { flex: 1 1 140px; min-width: 140px; display: grid; gap: 14px; align-content: start; }
.scale { display: grid; gap: 4px; color: var(--ink2); font-size: 12px; }
.bar { height: 8px; border-radius: 4px; }
.ends { display: flex; justify-content: space-between; color: var(--muted); font-size: 11px; font-variant-numeric: tabular-nums; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 18px; margin: 0 0 6px; }
.legend .row { display: inline-flex; align-items: baseline; gap: 6px; color: var(--ink2); }
.key { display: inline-block; width: 16px; height: 2px; border-radius: 1px; transform: translateY(-4px); }
.legend .val { color: var(--ink); font-weight: 600; font-variant-numeric: tabular-nums; }
#chart { display: block; width: 100%; touch-action: none; }
.tick { fill: var(--muted); font-size: 11px; font-variant-numeric: tabular-nums; }
.axis-title { fill: var(--ink2); font-size: 12px; }
.tip { position: absolute; pointer-events: none; background: var(--surface); border: 1px solid var(--border); border-radius: 8px; padding: 8px 10px; box-shadow: 0 4px 14px rgba(0,0,0,0.12); font-size: 12px; min-width: 150px; }
.tip .hd { color: var(--ink2); margin-bottom: 4px; font-variant-numeric: tabular-nums; }
.tip .r { display: grid; grid-template-columns: 16px auto 1fr; gap: 6px; align-items: center; }
.tip .r b { font-variant-numeric: tabular-nums; }
.tip .r span:last-child { color: var(--ink2); }
.note { color: var(--ink2); font-size: 12px; margin: 10px 0 0; max-width: 80ch; }
details { margin-top: 16px; }
summary { cursor: pointer; font-weight: 600; }
table { border-collapse: collapse; margin-top: 10px; font-size: 12px; font-variant-numeric: tabular-nums; }
th, td { padding: 3px 12px 3px 0; text-align: right; border-bottom: 1px solid var(--grid); }
th { color: var(--ink2); font-weight: 600; }
.tablewrap { overflow-x: auto; }
</style>
</head>
<body>
<main>
  <h1 id="title"></h1>
  <p class="sub" id="sub"></p>
  <div class="controls">
    <button id="play" type="button">Play</button>
    <div class="time">
      <input id="slider" type="range" min="0" step="0.5" value="0" aria-label="Time after the injection starts, in seconds">
      <span class="readout" id="readout">0.0 s</span>
    </div>
    <label>Speed <select id="speed"><option value="1">1×</option><option value="2">2×</option><option value="5" selected>5×</option><option value="10">10×</option></select></label>
    <label><input type="checkbox" id="show-artery" checked> Arteries</label>
    <label><input type="checkbox" id="show-vein" checked> Veins</label>
    <label><input type="checkbox" id="show-tissue" checked> Tissue</label>
  </div>
  <div class="panels">
    <section class="card legcard" aria-label="The leg's main vessels and tissue, coloured by iodine">
      <h2>Where the iodine is</h2>
      <div class="legwrap">
        <canvas id="leg" role="img" aria-label="Front view of the leg's main arteries, veins and tissue patches, coloured by their iodine at the time shown"></canvas>
        <div class="scales" id="scales"></div>
      </div>
      <p class="note" id="legnote"></p>
    </section>
    <section class="card chartcard">
      <h2>Iodine in the leg</h2>
      <div class="legend" id="legend"></div>
      <svg id="chart" tabindex="0" role="img" aria-label="Iodine in the leg's arteries, tissue and veins over time; arrow keys move the time"></svg>
      <div class="tip" id="tip" hidden></div>
      <p class="note" id="chartnote"></p>
    </section>
  </div>
  <details class="card">
    <summary>Data table</summary>
    <div class="tablewrap"><table id="table"></table></div>
  </details>
</main>
<script id="data" type="application/json">__DATA__</script>
<script>
(() => {
"use strict";
const D = JSON.parse(document.getElementById("data").textContent);
const b64 = s => { const bin = atob(s), u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return u; };
const VQ = b64(D.vein.q), TQ = b64(D.tis.q);
const NA = D.art.x1.length, NV = D.vein.x1.length, NT = D.tis.x.length;
const T = D.t, NF = T.length, FDT = T[1] - T[0], TEND = T[NF - 1];
const S = D.series, NS = S.keys.length;
const $ = id => document.getElementById(id);
let t = 0, playing = false, speed = 5, prevStamp = null;
const show = { artery: true, vein: true, tissue: true };

document.title = D.meta.title;
$("title").textContent = D.meta.title;
$("sub").textContent = D.meta.subtitle;
$("legnote").textContent = D.meta.legnote;
$("chartnote").textContent = D.meta.chartnote;

// ---------------------------------------------------------------- colours, from the theme
let C = {};
const rgb = hex => { const h = hex.trim().replace("#", ""); return [0, 2, 4].map(i => parseInt(h.slice(i, i + 2), 16)); };
const mix = (a, b, v) => { v = Math.max(0, Math.min(1, v)); return "rgb(" + Math.round(a[0] + (b[0] - a[0]) * v) + "," + Math.round(a[1] + (b[1] - a[1]) * v) + "," + Math.round(a[2] + (b[2] - a[2]) * v) + ")"; };
function readTheme() {
  const cs = getComputedStyle(document.documentElement), g = n => cs.getPropertyValue(n).trim();
  C = { surface: g("--surface"), leg: g("--leg"), empty: rgb(g("--empty")), artery: rgb(g("--artery")), vein: rgb(g("--vein")),
        tissue: rgb(g("--tissue")), ink: g("--ink"), ink2: g("--ink2"), muted: g("--muted"), grid: g("--grid"), axis: g("--axis"),
        series: S.keys.map(k => g("--" + k)) };
}
// magnitude: one hue from the empty grey, on a square-root scale so faint iodine still shows
const shade = (hue, v) => mix(C.empty, hue, Math.sqrt(Math.max(0, v)));

// ---------------------------------------------------------------- values at time t
function aifAt(s) {
  if (s <= 0) return 0;
  const a = D.aif, x = s / D.aif_dt, i = Math.floor(x);
  if (i >= a.length - 1) return a[a.length - 1];
  return a[i] + (x - i) * (a[i + 1] - a[i]);
}
function frame(tt) { const x = Math.min(Math.max(tt, 0) / FDT, NF - 1), f0 = Math.floor(x); return [f0, Math.min(f0 + 1, NF - 1), x - f0]; }
function qAt(Q, n, j, fr) { return ((1 - fr[2]) * Q[fr[0] * n + j] + fr[2] * Q[fr[1] * n + j]) / 255; }
function seriesAt(k, fr) { const v = S.values[k]; return (1 - fr[2]) * v[fr[0]] + fr[2] * v[fr[1]]; }
const fmtMg = v => v >= 100 ? Math.round(v).toLocaleString("en-US") : v >= 10 ? v.toFixed(1) : v.toFixed(2);

// ---------------------------------------------------------------- the leg
const canvas = $("leg"), ctx = canvas.getContext("2d"), B = D.bounds;
let scale = 1, cw = 0, ch = 0;
const X = x => (x - B.x0) * scale, Z = z => (B.z1 - z) * scale;
const dotR = D.tis.ml.map(ml => Math.max(1.4, Math.min(4.5, 1.15 * Math.cbrt(ml))));
function layoutLeg() {
  const card = canvas.closest(".card"), inner = card.clientWidth - 34, aspect = (B.z1 - B.z0) / (B.x1 - B.x0);
  let h = Math.max(360, window.innerHeight - 240), w = h / aspect;
  if (w > inner) { w = inner; h = w * aspect; }
  const dpr = window.devicePixelRatio || 1;
  canvas.style.width = w + "px"; canvas.style.height = h + "px";
  canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  scale = w / (B.x1 - B.x0); cw = w; ch = h;
}
function vessel(s, j, colour) {
  ctx.strokeStyle = colour;
  ctx.lineWidth = Math.max(0.6, 2 * s.r[j] * scale);
  ctx.beginPath(); ctx.moveTo(X(s.x1[j]), Z(s.z1[j])); ctx.lineTo(X(s.x2[j]), Z(s.z2[j])); ctx.stroke();
}
function drawLeg() {
  const fr = frame(t);
  ctx.fillStyle = C.surface; ctx.fillRect(0, 0, cw, ch);
  ctx.fillStyle = C.leg;
  const h = D.sil.h;
  for (const row of D.sil.rows) {
    for (let i = 1; i < row.length; i += 2) ctx.fillRect(X(row[i]), Z(row[0] + h), (row[i + 1] - row[i] + h) * scale + 0.5, h * scale + 0.5);
  }
  if (show.tissue) {
    for (let j = 0; j < NT; j++) {
      ctx.fillStyle = shade(C.tissue, qAt(TQ, NT, j, fr));
      ctx.beginPath(); ctx.arc(X(D.tis.x[j]), Z(D.tis.z[j]), dotR[j], 0, 2 * Math.PI); ctx.fill();
    }
  }
  ctx.lineCap = "round";
  if (show.vein) for (let j = 0; j < NV; j++) vessel(D.vein, j, shade(C.vein, qAt(VQ, NV, j, fr)));
  if (show.artery) for (let j = 0; j < NA; j++) vessel(D.art, j, shade(C.artery, aifAt(t - D.art.tau[j]) / D.aif_max));
}
function buildScales() {
  const box = $("scales"); box.textContent = "";
  const rows = [["artery", "Arteries, blood", D.aif_max, "mg/mL"], ["vein", "Veins, blood", D.vein.max, "mg/mL"], ["tissue", "Tissue patches", D.tis.max, "mg/mL"]];
  for (const [key, name, max, unit] of rows) {
    const stops = [0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1].map(p => shade(C[key], p) + " " + (100 * p) + "%").join(", ");
    const row = document.createElement("div"); row.className = "scale";
    const label = document.createElement("span"); label.textContent = name;
    const bar = document.createElement("div"); bar.className = "bar"; bar.style.background = "linear-gradient(to right, " + stops + ")";
    const ends = document.createElement("div"); ends.className = "ends";
    const a = document.createElement("span"); a.textContent = "0";
    const b = document.createElement("span"); b.textContent = (max >= 1 ? max.toFixed(1) : max.toFixed(3)) + " " + unit;
    ends.append(a, b); row.append(label, bar, ends); box.append(row);
  }
  const key = document.createElement("p"); key.className = "note"; key.style.margin = "0";
  key.textContent = "Grey: no iodine yet. Square-root scale, so faint iodine still shows.";
  box.append(key);
}

// ---------------------------------------------------------------- the chart
const svg = $("chart"), NSVG = "http://www.w3.org/2000/svg", tip = $("tip");
const M = { l: 58, r: 18, t: 14, b: 44 };
let W = 600, H = 360, ymax = 1, xs = v => v, ys = v => v, clip, now, hair, dots = [], overlay;
function el(name, attrs, parent) { const e = document.createElementNS(NSVG, name); for (const k in attrs) e.setAttribute(k, attrs[k]); (parent || svg).appendChild(e); return e; }
function niceTicks(max) {
  const raw = max / 4, p = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * p).find(s => s >= raw);
  const ticks = []; for (let v = 0; v <= Math.ceil(max / step) * step + 1e-9; v += step) ticks.push(v);
  return { top: ticks[ticks.length - 1], ticks };
}
function buildChart() {
  W = Math.max(280, svg.parentElement.clientWidth - 32); H = Math.round(Math.max(260, Math.min(430, W * 0.62)));
  svg.setAttribute("viewBox", "0 0 " + W + " " + H); svg.setAttribute("height", H);
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const nt = niceTicks(Math.max(...S.values.map(v => Math.max(...v)))); ymax = nt.top;
  xs = s => M.l + (W - M.l - M.r) * s / TEND; ys = v => H - M.b - (H - M.t - M.b) * v / ymax;
  for (const v of nt.ticks) {
    el("line", { x1: M.l, x2: W - M.r, y1: ys(v), y2: ys(v), stroke: v === 0 ? C.axis : C.grid, "stroke-width": 1 });
    const tx = el("text", { x: M.l - 8, y: ys(v) + 4, "text-anchor": "end", class: "tick" }); tx.textContent = v.toLocaleString("en-US");
  }
  for (let s = 0; s <= TEND + 1e-9; s += 20) { const tx = el("text", { x: xs(s), y: H - M.b + 18, "text-anchor": "middle", class: "tick" }); tx.textContent = s; }
  const xt = el("text", { x: (M.l + W - M.r) / 2, y: H - 8, "text-anchor": "middle", class: "axis-title" }); xt.textContent = "Time after the injection starts (s)";
  const yt = el("text", { x: 14, y: (M.t + H - M.b) / 2, "text-anchor": "middle", class: "axis-title", transform: "rotate(-90 14 " + (M.t + H - M.b) / 2 + ")" }); yt.textContent = "Iodine (mg)";
  const defs = el("defs", {}); const cp = el("clipPath", { id: "past" }, defs);
  clip = el("rect", { x: 0, y: 0, width: 0, height: H }, cp);
  for (let k = 0; k < NS; k++) {
    let d = ""; const v = S.values[k];
    for (let f = 0; f < NF; f++) d += (f ? "L" : "M") + xs(T[f]).toFixed(1) + "," + ys(v[f]).toFixed(1);
    el("path", { d, fill: "none", stroke: C.series[k], "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round", opacity: 0.22 });
    el("path", { d, fill: "none", stroke: C.series[k], "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round", "clip-path": "url(#past)" });
  }
  now = el("line", { y1: M.t, y2: H - M.b, stroke: C.ink2, "stroke-width": 1 });
  hair = el("line", { y1: M.t, y2: H - M.b, stroke: C.muted, "stroke-width": 1, visibility: "hidden" });
  dots = S.keys.map((_, k) => el("circle", { r: 4, fill: C.series[k], stroke: C.surface, "stroke-width": 2 }));
  overlay = el("rect", { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, fill: "transparent", style: "cursor: crosshair" });
  overlay.addEventListener("pointermove", hover);
  overlay.addEventListener("pointerleave", () => { tip.hidden = true; hair.setAttribute("visibility", "hidden"); });
  overlay.addEventListener("pointerdown", e => { pause(); t = timeAt(e); render(); hover(e); });
}
function timeAt(e) {
  const r = svg.getBoundingClientRect(), x = (e.clientX - r.left) * W / r.width;
  const s = (x - M.l) / (W - M.l - M.r) * TEND;
  return Math.max(0, Math.min(TEND, Math.round(s / FDT) * FDT));
}
function hover(e) {
  const s = timeAt(e), fr = frame(s), x = xs(s);
  hair.setAttribute("x1", x); hair.setAttribute("x2", x); hair.setAttribute("visibility", "visible");
  tip.textContent = "";
  const hd = document.createElement("div"); hd.className = "hd"; hd.textContent = s.toFixed(1) + " s"; tip.append(hd);
  for (let k = 0; k < NS; k++) {
    const r = document.createElement("div"); r.className = "r";
    const key = document.createElement("span"); key.className = "key"; key.style.background = C.series[k]; key.style.transform = "none";
    const b = document.createElement("b"); b.textContent = fmtMg(seriesAt(k, fr)) + " mg";
    const n = document.createElement("span"); n.textContent = S.names[k];
    r.append(key, b, n); tip.append(r);
  }
  tip.hidden = false;
  const card = svg.parentElement.getBoundingClientRect(), sr = svg.getBoundingClientRect();
  let left = e.clientX - card.left + 14, top = e.clientY - card.top + 14;
  if (left + tip.offsetWidth > card.width - 8) left = e.clientX - card.left - tip.offsetWidth - 14;
  if (top + tip.offsetHeight > sr.bottom - card.top) top = sr.bottom - card.top - tip.offsetHeight - 4;
  tip.style.left = left + "px"; tip.style.top = top + "px";
}
function buildLegend() {
  const box = $("legend"); box.textContent = "";
  S.keys.forEach((_, k) => {
    const row = document.createElement("span"); row.className = "row";
    const key = document.createElement("span"); key.className = "key"; key.style.background = C.series[k];
    const n = document.createElement("span"); n.textContent = S.names[k];
    const v = document.createElement("span"); v.className = "val"; v.id = "val" + k;
    row.append(key, n, v); box.append(row);
  });
}
function updateChart() {
  const fr = frame(t), x = xs(t);
  clip.setAttribute("width", Math.max(0, x));
  now.setAttribute("x1", x); now.setAttribute("x2", x);
  for (let k = 0; k < NS; k++) {
    const v = seriesAt(k, fr);
    dots[k].setAttribute("cx", x); dots[k].setAttribute("cy", ys(v));
    $("val" + k).textContent = fmtMg(v) + " mg";
  }
}
function buildTable() {
  const tb = $("table"); tb.textContent = "";
  const head = tb.createTHead().insertRow();
  for (const h of ["Time (s)"].concat(S.names.map(n => n + " (mg)"))) { const th = document.createElement("th"); th.textContent = h; head.append(th); }
  const body = tb.createTBody(), stride = Math.max(1, Math.round(5 / FDT));
  for (let f = 0; f < NF; f += stride) {
    const r = body.insertRow(); r.insertCell().textContent = T[f].toFixed(0);
    for (let k = 0; k < NS; k++) r.insertCell().textContent = fmtMg(S.values[k][f]);
  }
}

// ---------------------------------------------------------------- time, controls
const slider = $("slider"), play = $("play");
slider.max = TEND;
function render() {
  slider.value = t; $("readout").textContent = t.toFixed(1) + " s";
  drawLeg(); updateChart();
}
function pause() { playing = false; play.textContent = "Play"; }
function tick(stamp) {
  if (!playing) { prevStamp = null; return; }
  if (prevStamp !== null) t = Math.min(TEND, t + (stamp - prevStamp) / 1000 * speed);
  prevStamp = stamp; render();
  if (t >= TEND) { pause(); return; }
  requestAnimationFrame(tick);
}
play.addEventListener("click", () => {
  if (playing) { pause(); return; }
  if (t >= TEND) t = 0;
  playing = true; play.textContent = "Pause"; prevStamp = null; requestAnimationFrame(tick);
});
slider.addEventListener("input", () => { pause(); t = parseFloat(slider.value); render(); });
$("speed").addEventListener("change", e => { speed = parseFloat(e.target.value); });
for (const k of ["artery", "vein", "tissue"]) $("show-" + k).addEventListener("change", e => { show[k] = e.target.checked; drawLeg(); });
svg.addEventListener("keydown", e => {
  const step = { ArrowRight: 1, ArrowLeft: -1, PageUp: 10, PageDown: -10 }[e.key];
  if (step !== undefined) { pause(); t = Math.max(0, Math.min(TEND, t + step)); render(); e.preventDefault(); }
  else if (e.key === "Home" || e.key === "End") { pause(); t = e.key === "Home" ? 0 : TEND; render(); e.preventDefault(); }
});
function rebuild() { readTheme(); layoutLeg(); buildScales(); buildLegend(); buildChart(); render(); }
let resizeTimer = null;
window.addEventListener("resize", () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(rebuild, 120); });
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", rebuild);
new MutationObserver(rebuild).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
buildTable(); rebuild();
})();
</script>
</body>
</html>
"""

"Writes the page with data `D` to `path`."
function write_view(path, D)
    html = replace(VIEW_TEMPLATE, "__DATA__" => view_json(D))
    write(path, html)
    return path
end
