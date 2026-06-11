/* ============================================================
   Quantum Lab — rendering, UI, live diagnostics
   (qm-engine.js = TDSE solver, qm-scenarios.js = experiments)
   ============================================================ */

"use strict";

const $ = id => document.getElementById(id);

let S = null;               // current QMState
let SC = null;              // current scenario object
let SCP = {};               // current parameter values
let paused = false;
let E0 = null, N0 = null;
const consHist = [];        // {dN, dE}
const uncHist = [];         // Δx·Δp over time
const HIST_MAX = 500;
let eigenResult = null;

// ---------------- canvas ----------------

const cv = $("sim");
const ctx = cv.getContext("2d");
let W = 0, H = 0, DPR = 1;
let img = null;             // ImageData for 2D fields

function resize() {
  DPR = Math.min(2, window.devicePixelRatio || 1);
  W = window.innerWidth; H = window.innerHeight;
  cv.width = W * DPR; cv.height = H * DPR;
  cv.style.width = W + "px"; cv.style.height = H + "px";
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
}
window.addEventListener("resize", resize);
resize();

const VIEW_W = () => W - 380;

// phase → hue → rgb (fast HSV with s=0.85)
function phaseRGB(phase, v) {
  const h = (phase / (2 * Math.PI) + 0.5) * 6;
  const i = Math.floor(h) % 6;
  const f = h - Math.floor(h);
  const p = v * 0.15, q = v * (1 - 0.85 * f), t = v * (1 - 0.85 * (1 - f));
  switch (i) {
    case 0: return [v, t, p];
    case 1: return [q, v, p];
    case 2: return [p, v, t];
    case 3: return [p, q, v];
    case 4: return [t, p, v];
    default: return [v, p, q];
  }
}

// ---------------- 1D rendering ----------------

let peak1d = 1e-9;

function draw1D() {
  const w = VIEW_W(), h = H;
  const pad = 60, plotH = h - 2 * pad;
  const n = S.n;
  // half-domain scenarios (mirror-trick wall): show only the physical x ≥ 0
  const i0 = SC.halfDomain ? n / 2 : 0;
  const cnt = n - i0;
  const x0World = S.x(i0);
  const spanWorld = S.L - (x0World - S.x(0));

  // axes / potential scaling
  let vMin = Infinity, vMax = -Infinity;
  for (let i = 0; i < n; i++) { vMin = Math.min(vMin, S.V[i]); vMax = Math.max(vMax, S.V[i]); }
  const vSpanV = Math.max(vMax - vMin, 1e-9);

  // |ψ|² with auto-normalized peak (smoothed so it doesn't pulse)
  let pk = 1e-12;
  for (let i = 0; i < n; i++) pk = Math.max(pk, S.re[i] ** 2 + S.im[i] ** 2);
  peak1d = Math.max(pk, peak1d * 0.985);

  const X = i => (i - i0) / cnt * w;
  const Ypsi = d => h - pad - (d / peak1d) * plotH * 0.82;
  const Yv = v => h - pad - ((v - vMin) / vSpanV) * plotH * 0.55;

  // potential
  ctx.strokeStyle = "rgba(124,155,255,.5)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  for (let i = i0; i < n; i++) {
    const y = Yv(Math.min(S.V[i], vMin + vSpanV));
    i > i0 ? ctx.lineTo(X(i), y) : ctx.moveTo(X(i), y);
  }
  ctx.stroke();

  // |ψ|² as phase-colored columns
  const cols = Math.min(w, 900);
  for (let cx = 0; cx < cols; cx++) {
    const i = Math.min(n - 1, i0 + Math.round(cx / cols * cnt));
    const d = S.re[i] ** 2 + S.im[i] ** 2;
    if (d / peak1d < 1e-4) continue;
    const ph = Math.atan2(S.im[i], S.re[i]);
    const [r, g, b] = phaseRGB(ph, 1);
    const x = cx / cols * w, y = Ypsi(d);
    ctx.fillStyle = `rgba(${(r * 255) | 0},${(g * 255) | 0},${(b * 255) | 0},0.75)`;
    ctx.fillRect(x, y, w / cols + 0.5, h - pad - y);
  }
  // envelope
  ctx.strokeStyle = "#fff";
  ctx.lineWidth = 1.4;
  ctx.beginPath();
  for (let i = i0; i < n; i++) {
    const y = Ypsi(S.re[i] ** 2 + S.im[i] ** 2);
    i > i0 ? ctx.lineTo(X(i), y) : ctx.moveTo(X(i), y);
  }
  ctx.stroke();

  // analytic overlay: free-packet envelope
  if (SC === QM_SCENARIOS.free) {
    const m = S.moments();
    const sp = SC.sigmaPred(SCP, S.t);
    const A = peak1d;                          // align peaks for shape comparison
    ctx.strokeStyle = "#ffd166";
    ctx.setLineDash([5, 5]);
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const x = S.x(i);
      const d = A * Math.exp(-((x - m.x) ** 2) / (2 * sp * sp));
      const y = Ypsi(d);
      i ? ctx.lineTo(X(i), y) : ctx.moveTo(X(i), y);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // ⟨x⟩ marker
  const m = S.moments(SC.halfDomain ? 0 : -Infinity);
  const mx = (m.x - x0World) / spanWorld * w;
  ctx.strokeStyle = "rgba(255,255,255,.6)";
  ctx.setLineDash([3, 4]);
  ctx.beginPath(); ctx.moveTo(mx, pad); ctx.lineTo(mx, h - pad); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "rgba(255,255,255,.7)";
  ctx.font = "11px ui-monospace, monospace";
  ctx.fillText("⟨x⟩", mx + 4, pad + 12);

  // x axis
  ctx.strokeStyle = "rgba(255,255,255,.25)";
  ctx.beginPath(); ctx.moveTo(0, h - pad); ctx.lineTo(w, h - pad); ctx.stroke();
  ctx.fillStyle = "#8890b8";
  const tickStep = spanWorld > 100 ? 20 : 5;
  for (let xt = Math.ceil(x0World / tickStep) * tickStep; xt <= S.L / 2; xt += tickStep) {
    const px = (xt - x0World) / spanWorld * w;
    ctx.fillText(String(xt), px - 8, h - pad + 16);
  }
}

// ---------------- 2D rendering ----------------

let peak2d = 1e-9;

function draw2D() {
  const { nx, ny } = S;
  if (!img || img.width !== nx) img = ctx.createImageData(nx, ny);
  const data = img.data;

  let pk = 1e-12;
  for (let i = 0; i < S.size; i++) pk = Math.max(pk, S.re[i] ** 2 + S.im[i] ** 2);
  peak2d = Math.max(pk, peak2d * 0.97);

  for (let i = 0; i < S.size; i++) {
    const d = S.re[i] ** 2 + S.im[i] ** 2;
    const v = Math.min(1, Math.pow(d / peak2d, 0.42));
    const [r, g, b] = v > 0.004 ? phaseRGB(Math.atan2(S.im[i], S.re[i]), v) : [0, 0, 0];
    const o = i * 4;
    // potential shown as dim blue-gray
    const vp = S.V[i];
    const wall = vp > 1 ? 0.35 : Math.max(0, -vp) * 0.04;
    data[o] = Math.min(255, ((r * 0.97 + wall * 0.25) * 255) | 0);
    data[o + 1] = Math.min(255, ((g * 0.97 + wall * 0.3) * 255) | 0);
    data[o + 2] = Math.min(255, ((b * 0.97 + wall * 0.55) * 255) | 0);
    data[o + 3] = 255;
  }

  const w = VIEW_W(), h = H;
  const side = Math.min(w - 40, h - 80);
  const ox = (w - side) / 2, oy = (h - side) / 2;

  // draw via offscreen canvas scale-up
  if (!draw2D.oc || draw2D.oc.width !== nx) {
    draw2D.oc = document.createElement("canvas");
    draw2D.oc.width = nx; draw2D.oc.height = ny;
  }
  draw2D.oc.getContext("2d").putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(draw2D.oc, ox, oy, side, side);
  ctx.strokeStyle = "rgba(124,155,255,.35)";
  ctx.strokeRect(ox, oy, side, side);

  const toPx = (x, y) => [
    ox + (x + S.Lx / 2) / S.Lx * side,
    oy + (y + S.Ly / 2) / S.Ly * side,
  ];

  // classical overlay (Ehrenfest)
  if (SC._cl) {
    const c = SC._cl;
    ctx.strokeStyle = "rgba(255,255,255,.8)";
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let k = 0; k < c.trail.length / 2; k++) {
      const [px, py] = toPx(c.trail[k * 2], c.trail[k * 2 + 1]);
      k ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
    }
    ctx.stroke();
    const [px, py] = toPx(c.x, c.y);
    ctx.fillStyle = "#fff";
    ctx.beginPath(); ctx.arc(px, py, 4, 0, 7); ctx.fill();
    ctx.font = "11px ui-monospace, monospace";
    ctx.fillText("classical", px + 8, py - 6);
    const m = S.moments();
    const [qx, qy] = toPx(m.x, m.y);
    ctx.strokeStyle = "#ffd166";
    ctx.beginPath(); ctx.arc(qx, qy, 6, 0, 7); ctx.stroke();
    ctx.fillStyle = "#ffd166";
    ctx.fillText("⟨r⟩", qx + 9, qy + 4);
  }

  // screen strip (double slit)
  if (SC._screen) {
    const sx = toPx(SC.screenX, 0)[0];
    ctx.strokeStyle = "rgba(255,209,102,.5)";
    ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(sx, oy); ctx.lineTo(sx, oy + side); ctx.stroke();
    ctx.setLineDash([]);
    let mx = 1e-12;
    for (const v of SC._screen) mx = Math.max(mx, v);
    ctx.strokeStyle = "#ffd166";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (let j = 0; j < S.ny; j++) {
      const px = ox + side + 6 + (SC._screen[j] / mx) * 60;
      const py = oy + j / S.ny * side;
      j ? ctx.lineTo(px, py) : ctx.moveTo(px, py);
    }
    ctx.stroke();
  }

  // axis labels
  ctx.fillStyle = "#8890b8";
  ctx.font = "11px ui-monospace, monospace";
  const ax = SC.axes || ["x", "y"];
  ctx.fillText(ax[0], ox + side / 2 - 20, oy + side + 16);
  ctx.save();
  ctx.translate(ox - 8, oy + side / 2 + 30);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(ax[1], 0, 0);
  ctx.restore();
}

// ---------------- graphs ----------------

const gCons = $("gCons").getContext("2d");
const gUnc = $("gUnc").getContext("2d");
const gScen = $("gScen").getContext("2d");

function drawConsGraph() {
  const w = 320, h = 90;
  gCons.clearRect(0, 0, w, h);
  if (consHist.length < 2) return;
  const X = i => i / (HIST_MAX - 1) * w;
  const Y = v => {
    const lg = Math.log10(Math.abs(v) + 1e-16);
    return Math.max(2, Math.min(h - 2, h - (lg + 16) / 16 * h));
  };
  for (const [key, col] of [["dN", "#5ce8c5"], ["dE", "#ffd166"]]) {
    gCons.strokeStyle = col; gCons.lineWidth = 1.3;
    gCons.beginPath();
    consHist.forEach((s, i) => i ? gCons.lineTo(X(i), Y(s[key])) : gCons.moveTo(X(i), Y(s[key])));
    gCons.stroke();
  }
  gCons.fillStyle = "rgba(255,255,255,.4)";
  gCons.font = "9px ui-monospace, monospace";
  for (const lg of [-4, -8, -12]) {
    const y = h - (lg + 16) / 16 * h;
    gCons.fillRect(0, y, w, 0.5);
    gCons.fillText("1e" + lg, 3, y - 2);
  }
}

function drawUncGraph() {
  const w = 320, h = 90;
  gUnc.clearRect(0, 0, w, h);
  if (uncHist.length < 2) return;
  let hi = 0.6;
  for (const v of uncHist) hi = Math.max(hi, v);
  const X = i => i / (HIST_MAX - 1) * w;
  const Y = v => h - 4 - (v / hi) * (h - 12);
  gUnc.strokeStyle = "rgba(255,255,255,.4)";
  gUnc.setLineDash([4, 4]);
  gUnc.beginPath(); gUnc.moveTo(0, Y(0.5)); gUnc.lineTo(w, Y(0.5)); gUnc.stroke();
  gUnc.setLineDash([]);
  gUnc.fillStyle = "rgba(255,255,255,.5)";
  gUnc.font = "9px ui-monospace, monospace";
  gUnc.fillText("ħ/2", 4, Y(0.5) - 3);
  gUnc.strokeStyle = "#c792ea"; gUnc.lineWidth = 1.5;
  gUnc.beginPath();
  uncHist.forEach((v, i) => i ? gUnc.lineTo(X(i), Y(v)) : gUnc.moveTo(X(i), Y(v)));
  gUnc.stroke();
}

function drawScenGraph() {
  const w = 320, h = 110;
  gScen.clearRect(0, 0, w, h);
  const cfg = SC.graph;
  if (!cfg) return;
  gScen.fillStyle = "rgba(255,255,255,.5)";
  gScen.font = "9px ui-monospace, monospace";
  gScen.fillText(cfg.title, 4, 10);

  if (cfg.screen && SC._screen) {
    let mx = 1e-12;
    for (const v of SC._screen) mx = Math.max(mx, v);
    gScen.strokeStyle = "#ffd166"; gScen.lineWidth = 1.4;
    gScen.beginPath();
    SC._screen.forEach((v, j) => {
      const x = j / SC._screen.length * w, y = h - 4 - (v / mx) * (h - 20);
      j ? gScen.lineTo(x, y) : gScen.moveTo(x, y);
    });
    gScen.stroke();
    return;
  }
  const hist = SC._hist;
  if (!hist || hist.length < 2) return;
  let lo = Infinity, hi = -Infinity;
  for (const r of hist) {
    for (let k = 1; k <= 2; k++) if (r[k] !== null && r[k] !== undefined) {
      lo = Math.min(lo, r[k]); hi = Math.max(hi, r[k]);
    }
  }
  const pad = (hi - lo) * 0.1 + 1e-12; lo -= pad; hi += pad;
  const X = i => i / (hist.length - 1) * w;
  const Y = v => h - 4 - (v - lo) / (hi - lo) * (h - 20);
  // series 2 = analytic (dashed gold), series 1 = measured (solid)
  if (cfg.series >= 2) {
    gScen.strokeStyle = "#ffd166"; gScen.setLineDash([4, 4]); gScen.lineWidth = 1.2;
    gScen.beginPath();
    hist.forEach((r, i) => { if (r[2] !== null) i ? gScen.lineTo(X(i), Y(r[2])) : gScen.moveTo(X(i), Y(r[2])); });
    gScen.stroke();
    gScen.setLineDash([]);
  }
  gScen.strokeStyle = "#5ce8c5"; gScen.lineWidth = 1.5;
  gScen.beginPath();
  hist.forEach((r, i) => i ? gScen.lineTo(X(i), Y(r[1])) : gScen.moveTo(X(i), Y(r[1])));
  gScen.stroke();
}

// ---------------- diagnostics ----------------

let frameCount = 0;

function updateDiagnostics() {
  const norm = S.norm();
  const en = S.energy();
  const m = S.moments(SC.halfDomain ? 0 : -Infinity);
  const sp = S.spectral();
  if (N0 === null) { N0 = norm; E0 = en.E; }
  consHist.push({ dN: (norm - N0) / N0, dE: Math.abs(E0) > 1e-12 ? (en.E - E0) / Math.abs(E0) : 0 });
  if (consHist.length > HIST_MAX) consHist.shift();
  const unc = Math.sqrt(m.varX) * Math.sqrt(sp.varPx);
  uncHist.push(unc);
  if (uncHist.length > HIST_MAX) uncHist.shift();

  $("readout").textContent =
    `t       = ${S.t.toFixed(2)}\n` +
    `norm    = ${norm.toFixed(10)}${S.mask ? "  (absorbers on)" : ""}\n` +
    `⟨H⟩     = ${en.E.toFixed(6)}   ⟨T⟩=${en.T.toFixed(3)} ⟨V⟩=${en.V.toFixed(3)}\n` +
    `⟨x⟩,⟨p⟩ = ${m.x.toFixed(3)}, ${sp.px.toFixed(3)}\n` +
    `Δx·Δp   = ${unc.toFixed(4)}   (≥ ħ/2 = 0.5)`;

  const lines = SC.measure ? SC.measure(S, SCP).filter(Boolean) : [];
  $("scenMeasure").textContent = lines.join("\n");

  // eigen table
  if (eigenResult && SC.eigen) {
    let t = "state   E (sim)      E exact     " + (SC.eigen.real ? "real units" : "") + "\n";
    eigenResult.energies.forEach((E, n) => {
      const ex = SC.eigen.exact(n, SCP);
      t += `n=${n + 1}     ${E.toFixed(5)}     ${ex.toFixed(5)}     ${SC.eigen.real ? SC.eigen.real(E) : ""}\n`;
    });
    $("eigenTable").textContent = t;
    $("eigenBox").style.display = "";
  } else {
    $("eigenBox").style.display = "none";
  }
}

// ---------------- scenario loading ----------------

function loadScenario(name) {
  SC = QM_SCENARIOS[name];
  SCP = {};
  for (const p of SC.params || []) SCP[p.key] = p.val;
  buildState();
  buildParamUI();
  $("scenInfo").textContent = SC.info;
  $("eqV").textContent = "";
}

function buildState() {
  S = new QMState(SC.dim === 1
    ? { dim: 1, n: SC.n, L: SC.L }
    : { dim: 2, nx: SC.nx, ny: SC.ny, Lx: SC.Lx, Ly: SC.Ly });
  SC.build(S, SCP);
  E0 = N0 = null;
  consHist.length = 0; uncHist.length = 0;
  peak1d = peak2d = 1e-9;
  eigenResult = null;
  if (SC.eigen) {
    // synchronous relaxation; fast at these grid sizes
    eigenResult = solveEigen(S, SC.eigen.count, { dtau: 2e-3 });
  }
  $("sDt").value = Math.log10(SC.dt);
  $("sPF").value = SC.perFrame;
  updateOutputs();
}

function buildParamUI() {
  const box = $("params");
  box.innerHTML = "";
  for (const p of SC.params || []) {
    const lab = document.createElement("label");
    lab.innerHTML = `${p.label} <output>${p.val}</output>`;
    const inp = document.createElement("input");
    inp.type = "range"; inp.min = p.min; inp.max = p.max; inp.step = p.step; inp.value = p.val;
    inp.addEventListener("input", () => {
      SCP[p.key] = parseFloat(inp.value);
      lab.querySelector("output").textContent = inp.value;
      buildState();           // parameter changes restart the experiment
    });
    box.appendChild(lab);
    box.appendChild(inp);
  }
}

// ---------------- controls ----------------

let dt = 0.02, perFrame = 4;

function updateOutputs() {
  dt = Math.pow(10, parseFloat($("sDt").value));
  perFrame = Math.round(parseFloat($("sPF").value));
  $("oDt").textContent = dt.toExponential(1);
  $("oPF").textContent = perFrame;
}
$("sDt").addEventListener("input", updateOutputs);
$("sPF").addEventListener("input", updateOutputs);

$("preset").onchange = e => loadScenario(e.target.value);
const btnPlay = $("btnPlay");
function setPaused(p) { paused = p; btnPlay.textContent = paused ? "▶ Play" : "⏸ Pause"; }
btnPlay.onclick = () => setPaused(!paused);
$("btnStep").onclick = () => { setPaused(true); stepFrame(); };
$("btnReset").onclick = () => buildState();
$("btnMeasure").onclick = () => {
  // projective position measurement: sample x* from |ψ|², collapse to a
  // narrow Gaussian around it (finite detector resolution)
  if (S.dim !== 1) return;
  const weights = new Float64Array(S.n);
  let tot = 0;
  for (let i = 0; i < S.n; i++) { weights[i] = S.re[i] ** 2 + S.im[i] ** 2; tot += weights[i]; }
  let r = Math.random() * tot, idx = 0;
  for (let i = 0; i < S.n; i++) { r -= weights[i]; if (r <= 0) { idx = i; break; } }
  const xStar = S.x(idx), w = 4 * S.dx * 6;
  for (let i = 0; i < S.n; i++) {
    const g = Math.exp(-((S.x(i) - xStar) ** 2) / (4 * w * w));
    S.re[i] *= g; S.im[i] *= g;
  }
  S.normalize();
  E0 = N0 = null;            // measurement injects energy — reset baselines
};

window.addEventListener("keydown", e => {
  if (e.code === "Space" && !["INPUT", "SELECT"].includes(e.target.tagName)) {
    e.preventDefault(); setPaused(!paused);
  }
});

// ---------------- main loop ----------------

function stepFrame() {
  S.step(dt, perFrame);
  if (SC.tick) SC.tick(S, SCP, dt * perFrame);
}

function frame() {
  if (!paused) stepFrame();
  ctx.fillStyle = "#05060e";
  ctx.fillRect(0, 0, W, H);
  if (S.dim === 1) draw1D(); else draw2D();

  if (frameCount++ % 3 === 0) {
    updateDiagnostics();
    drawConsGraph();
    drawUncGraph();
    drawScenGraph();
    $("eqPsi").innerHTML =
      `<b>iħ ∂ψ/∂t</b> = −(ħ²/2m)∇²ψ + V·ψ &nbsp;·&nbsp; split-step Fourier (unitary) · ` +
      `${S.dim === 1 ? `N=${S.n}` : `${S.nx}×${S.ny}`} grid · ħ = m = 1`;
  }
  requestAnimationFrame(frame);
}

loadScenario("bouncer");
requestAnimationFrame(frame);
