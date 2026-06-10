/* ============================================================
   Gravity Lab — rendering, UI, measurement displays
   (physics.js = simulation core, presets.js = scenarios)
   ============================================================ */

"use strict";

// ---------------- view state ----------------

const V = {
  sheet: false,
  warp: 1,
  grid: true,
  trails: true,
  glow: true,
  follow: false,
  rotFrame: false,          // co-rotating frame of the two heaviest bodies
  yaw: 0.6, pitch: 1.05,
  scale: 120,               // px per AU
  cx: 0, cy: 0,             // camera center (view coords)
  launch: false,
  newMass: 1e-6,            // M☉
};

let paused = false;
let selected = null;
let curFrame = null;        // rotating-frame data for this render frame

let E0 = null, L0 = null, Lscale = 1;
const energyHist = [];      // {ke, pe, e}
const consHist = [];        // {dE, dL, dP}  (relative drifts)
const HIST_MAX = 400;
const TRAIL_MAX = 500;
// trail samples are stride-6: x, y, Φ, frame θ, frame bary x, y — so trails
// can be replayed correctly in the rotating frame and on the 3D sheet

let recording = [];         // [{t, s:[[id,x,y,vx,vy],...]}]
let recStride = 1, recCount = 0;

const ghostRays = [];       // finished light-ray paths, fading out
onPhotonDied = trail => {
  ghostRays.push({ trail, alpha: 0.55 });
  if (ghostRays.length > 40) ghostRays.shift();
};

// ---------------- per-frame stepping ----------------

onBodyRemoved = (gone, into) => { if (selected === gone) selected = into; };

function resetBaselines() {
  E0 = null; L0 = null;
  energyHist.length = 0; consHist.length = 0;
}

// lightweight rotating-frame info (no Lagrange root-finding)
function frameInfo() {
  if (bodies.length < 2) return { theta: 0, bx: 0, by: 0 };
  let b1 = null, b2 = null;
  for (const b of bodies) {
    if (!b1 || b.m > b1.m) { b2 = b1; b1 = b; }
    else if (!b2 || b.m > b2.m) b2 = b;
  }
  const M = b1.m + b2.m;
  return {
    theta: Math.atan2(b2.y - b1.y, b2.x - b1.x),
    bx: (b1.x * b1.m + b2.x * b2.m) / M,
    by: (b1.y * b1.m + b2.y * b2.m) / M,
  };
}

function pushTrail(b, f) {
  b.trail.push(b.x, b.y, potentialAt(b.x, b.y, b), f.theta, f.bx, f.by);
  if (b.trail.length > TRAIL_MAX * 6) b.trail.splice(0, b.trail.length - TRAIL_MAX * 6);
  b._tx = b.x; b._ty = b.y;
}

function stepFrame() {
  // trails are sampled by distance moved (not per frame) so fast perihelion
  // passages stay smooth; the threshold tracks the current zoom level
  const thresh2 = (viewSpan() / 120) ** 2;
  for (let s = 0; s < P.sub; s++) {
    advance(P.dt);
    let f = null;
    for (const b of bodies) {
      const dx = b.x - (b._tx ?? 1e99), dy = b.y - (b._ty ?? 0);
      if (dx * dx + dy * dy > thresh2) pushTrail(b, f = f || frameInfo());
    }
  }
  if (P.merge) mergeCollisions();
  if (P.pw) captureByBH();

  if (photons.length) {
    const span = viewSpan();
    stepPhotons(P.dt * P.sub, (span * 4) ** 2, span / 220);
  }

  const en = energies();
  const L = angularMomentum();
  const [px, py] = momentum();
  if (E0 === null && bodies.length) {
    E0 = en.e; L0 = L;
    Lscale = Math.max(Math.abs(L), 1e-30);
    for (const b of bodies) Lscale = Math.max(Lscale, Math.abs(b.m) * Math.hypot(b.x, b.y) * Math.hypot(b.vx, b.vy));
  }
  energyHist.push(en);
  consHist.push({
    dE: Math.abs(E0) > 1e-300 ? (en.e - E0) / Math.abs(E0) : 0,
    dL: (L - (L0 || 0)) / Lscale,
    dP: Math.hypot(px, py) / (Lscale / Math.max(viewSpan(), 1e-9)),
  });
  if (energyHist.length > HIST_MAX) energyHist.shift();
  if (consHist.length > HIST_MAX) consHist.shift();

  // trajectory recording (auto-decimating ring)
  if (++recCount % recStride === 0) {
    recording.push({ t: time, s: bodies.map(b => [b.id, b.x, b.y, b.vx, b.vy]) });
    if (recording.length > 3000) {
      recording = recording.filter((_, i) => i % 2 === 0);
      recStride *= 2;
    }
  }
}

// ---------------- preset loading ----------------

function loadPreset(name) {
  // hard reset of physics + view modes (presets then override)
  bodies = []; photons = []; colorIdx = 0; time = 0; bodyIdSeq = 1;
  selected = null;
  keplerPoints.length = 0;

  ghostRays.length = 0;
  recording = []; recStride = 1; recCount = 0;
  resetBaselines();
  Object.assign(P, {
    Gmult: 1, cMult: 1, n: 2, eps: 0, dt: 0.5, sub: 12,
    integrator: "verlet", adaptive: true, pn1: false, pw: false, merge: false,
  });
  V.rotFrame = false; V.cx = 0; V.cy = 0; V.follow = false;

  const pr = PRESETS[name] || PRESETS.solar;
  pr.build();
  computeAccels();
  if (pr.select) selected = bodies.find(b => b.name === pr.select) || null;
  $("presetInfo").textContent = pr.info || "";
  syncUI();
}

// ---------------- canvas & projection ----------------

const cv = document.getElementById("sim");
const ctx = cv.getContext("2d");
let W = 0, H = 0, DPR = 1;

function resize() {
  DPR = Math.min(2, window.devicePixelRatio || 1);
  W = window.innerWidth; H = window.innerHeight;
  cv.width = W * DPR; cv.height = H * DPR;
  cv.style.width = W + "px"; cv.style.height = H + "px";
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
}
window.addEventListener("resize", resize);
resize();

const VIEW_W = () => W - 360;
const viewSpan = () => Math.max(VIEW_W(), 100) / V.scale * 0.55;

// inertial → rotating-frame view coords (identity when frame is off)
function viewXY(x, y, fh) {
  if (!V.rotFrame || !curFrame) return [x, y];
  const f = fh || curFrame.frame;
  const dx = x - f.bx, dy = y - f.by;
  const ct = Math.cos(-f.theta), st = Math.sin(-f.theta);
  return [curFrame.frame.bx + dx * ct - dy * st,
          curFrame.frame.by + dx * st + dy * ct];
}
// view coords → inertial
function frameInv(x, y) {
  if (!V.rotFrame || !curFrame) return [x, y];
  const f = curFrame.frame;
  const dx = x - f.bx, dy = y - f.by;
  const ct = Math.cos(f.theta), st = Math.sin(f.theta);
  return [f.bx + dx * ct - dy * st, f.by + dx * st + dy * ct];
}

// view coords → screen, always [sx, sy, depth]
function project(x, y, z) {
  const dx = x - V.cx, dy = y - V.cy;
  if (!V.sheet) return [VIEW_W() / 2 + dx * V.scale, H / 2 + dy * V.scale, 0];
  const cy_ = Math.cos(V.yaw), sy_ = Math.sin(V.yaw);
  const xr = dx * cy_ - dy * sy_;
  const yr = dx * sy_ + dy * cy_;
  const cp = Math.cos(V.pitch), sp = Math.sin(V.pitch);
  return [VIEW_W() / 2 + xr * V.scale, H / 2 + (yr * cp - (z || 0) * sp) * V.scale, yr];
}

function unproject(sx, sy) {     // screen → view coords on the z = 0 plane
  const px = (sx - VIEW_W() / 2) / V.scale;
  const py = (sy - H / 2) / V.scale;
  if (!V.sheet) return [px + V.cx, py + V.cy];
  const cp = Math.cos(V.pitch);
  const xr = px, yr = py / (cp || 1e-6);
  const cy_ = Math.cos(-V.yaw), sy_ = Math.sin(-V.yaw);
  return [xr * cy_ - yr * sy_ + V.cx, xr * sy_ + yr * cy_ + V.cy];
}

// potential used for grid color / sheet height at a *view* point
function gridPhiAt(vx, vy) {
  if (V.rotFrame && curFrame) {
    const [ix, iy] = frameInv(vx, vy);
    return effectivePotentialAt(ix, iy, curFrame.frame);
  }
  return potentialAt(vx, vy);
}

// self-normalizing sheet depth (visual only; physics is un-normalized)
let phiLo = -1, phiHi = 0, gridSpanCache = 1;
function depthT(phi) {
  return Math.max(0, Math.min(1, (phiHi - phi) / (phiHi - phiLo + 1e-300)));
}
function sheetZ(phi) {
  const maxDepth = Math.min(V.warp * gridSpanCache * 0.45,
                            H * 0.3 / (V.scale * Math.max(0.3, Math.sin(V.pitch))));
  return -maxDepth * Math.pow(depthT(phi), 0.75);
}

function bodyDrawRadius(b) {
  return Math.max(2.5, Math.min(16, b.radAU * V.scale * 15));
}

// ---------------- rendering ----------------

function heatColor(t) {
  t = Math.max(0, Math.min(1, t));
  const r = Math.round(40 + 215 * Math.pow(t, 1.4));
  const g = Math.round(60 + 130 * t * (1 - t) * 2.4);
  const b = Math.round(140 + 100 * (1 - t));
  return `rgb(${r},${g},${b})`;
}

function drawGrid() {
  const NG = V.sheet ? 36 : 30;
  const span = Math.max(VIEW_W(), H) / V.scale * (V.sheet ? 0.85 : 0.6);
  const step = (2 * span) / NG;
  const x0 = V.cx - span, y0 = V.cy - span;
  gridSpanCache = span;

  const phis = [];
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i <= NG; i++) {
    for (let j = 0; j <= NG; j++) {
      const phi = gridPhiAt(x0 + i * step, y0 + j * step);
      phis.push(phi);
      if (phi < lo) lo = phi;
      if (phi > hi) hi = phi;
    }
  }
  phiLo = lo; phiHi = hi;
  if (!V.grid) return;

  const bodyViews = bodies.map(b => viewXY(b.x, b.y));
  const pts = [];
  for (let i = 0; i <= NG; i++) {
    for (let j = 0; j <= NG; j++) {
      const k = i * (NG + 1) + j;
      const wx = x0 + i * step, wy = y0 + j * step;
      const t = depthT(phis[k]);
      let p;
      if (V.sheet) {
        p = project(wx, wy, sheetZ(phis[k]));
      } else {
        let ox = 0, oy = 0, dMin = Infinity;
        for (let bi = 0; bi < bodies.length; bi++) {
          const dx = bodyViews[bi][0] - wx, dy = bodyViews[bi][1] - wy;
          const d = Math.sqrt(dx * dx + dy * dy) + 1e-12;
          dMin = Math.min(dMin, d);
          const w = bodies[bi].m / Math.pow(d + P.eps, P.n - 1);
          ox += dx / d * w; oy += dy / d * w;
        }
        const len = Math.sqrt(ox * ox + oy * oy) + 1e-300;
        // capped at half the distance to the nearest body so nodes never overshoot
        const pull = Math.min(V.warp * step * 1.6 * Math.pow(t, 0.9), dMin * 0.5);
        p = project(wx + ox / len * pull, wy + oy / len * pull, 0);
      }
      p[3] = Math.sqrt(t);
      pts[k] = p;
    }
  }

  ctx.lineWidth = 1;
  const line = (a, b) => {
    const t = (a[3] + b[3]) / 2;
    ctx.strokeStyle = heatColor(t);
    ctx.globalAlpha = 0.16 + 0.5 * t;
    ctx.beginPath();
    ctx.moveTo(a[0], a[1]);
    ctx.lineTo(b[0], b[1]);
    ctx.stroke();
  };
  for (let i = 0; i <= NG; i++) {
    for (let j = 0; j <= NG; j++) {
      const k = i * (NG + 1) + j;
      if (i < NG) line(pts[k], pts[k + NG + 1]);
      if (j < NG) line(pts[k], pts[k + 1]);
    }
  }
  ctx.globalAlpha = 1;
}

function drawLagrangePoints() {
  if (!V.rotFrame || !curFrame) return;
  ctx.font = "11px ui-monospace, monospace";
  curFrame.points.forEach(([ix, iy], i) => {
    const [vx, vy] = viewXY(ix, iy);
    const p = V.sheet ? project(vx, vy, sheetZ(gridPhiAt(vx, vy))) : project(vx, vy, 0);
    ctx.strokeStyle = "#aef07a";
    ctx.globalAlpha = 0.9;
    ctx.beginPath(); ctx.arc(p[0], p[1], 5, 0, 7); ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(p[0] - 8, p[1]); ctx.lineTo(p[0] + 8, p[1]);
    ctx.moveTo(p[0], p[1] - 8); ctx.lineTo(p[0], p[1] + 8);
    ctx.stroke();
    ctx.fillStyle = "#aef07a";
    ctx.fillText("L" + (i + 1), p[0] + 8, p[1] - 8);
    ctx.globalAlpha = 1;
  });
}

function drawBHmarkers() {
  if (!P.pw) return;
  const bh = heaviest();
  if (!bh) return;
  const [vx, vy] = viewXY(bh.x, bh.y);
  const rs = schwarzschildR(bh.m);
  const rings = [[2 * rs, "rgba(255,107,138,.8)", "2 rₛ capture"],
                 [6 * Gv() * bh.m / (Cv() * Cv()), "rgba(255,209,102,.55)", "ISCO"]];
  ctx.setLineDash([5, 5]);
  for (const [r, col, label] of rings) {
    ctx.strokeStyle = col;
    ctx.beginPath();
    for (let a = 0; a <= 64; a++) {
      const x = vx + r * Math.cos(a / 64 * 2 * Math.PI);
      const y = vy + r * Math.sin(a / 64 * 2 * Math.PI);
      const p = V.sheet ? project(x, y, sheetZ(gridPhiAt(x, y))) : project(x, y, 0);
      a ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]);
    }
    ctx.stroke();
    const pl = project(vx + r * 0.71, vy - r * 0.71, 0);
    ctx.fillStyle = col;
    ctx.font = "10px ui-monospace, monospace";
    if (!V.sheet) ctx.fillText(label, pl[0] + 4, pl[1] - 4);
  }
  ctx.setLineDash([]);
}

function trailPoint(b, k) {
  const o = k * 6;
  const tx = b.trail[o], ty = b.trail[o + 1];
  let vx = tx, vy = ty;
  if (V.rotFrame && curFrame) {
    [vx, vy] = viewXY(tx, ty, { theta: b.trail[o + 3], bx: b.trail[o + 4], by: b.trail[o + 5] });
  }
  return V.sheet ? project(vx, vy, sheetZ(b.trail[o + 2])) : project(vx, vy, 0);
}

function drawBodies() {
  const order = [...bodies];
  if (V.sheet) order.sort((a, b) => {
    const va = viewXY(a.x, a.y), vb = viewXY(b.x, b.y);
    return project(va[0], va[1], 0)[2] - project(vb[0], vb[1], 0)[2];
  });

  for (const b of order) {
    if (V.trails && b.trail.length > 12) {
      const n = b.trail.length / 6;
      ctx.lineWidth = 1.5;
      ctx.strokeStyle = b.color;
      ctx.beginPath();
      for (let k = 0; k < n; k++) {
        const pp = trailPoint(b, k);
        k ? ctx.lineTo(pp[0], pp[1]) : ctx.moveTo(pp[0], pp[1]);
      }
      ctx.globalAlpha = 0.35;
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    const phi = potentialAt(b.x, b.y, b);
    const [vx, vy] = viewXY(b.x, b.y);
    const p = V.sheet ? project(vx, vy, sheetZ(V.rotFrame ? gridPhiAt(vx, vy) : phi))
                      : project(vx, vy, 0);
    const rad = bodyDrawRadius(b);

    if (V.glow) {
      const g = ctx.createRadialGradient(p[0], p[1], 0, p[0], p[1], rad * 5);
      g.addColorStop(0, b.color);
      g.addColorStop(0.25, b.color + "55");
      g.addColorStop(1, "transparent");
      ctx.fillStyle = g;
      ctx.beginPath(); ctx.arc(p[0], p[1], rad * 5, 0, 7); ctx.fill();
    }
    ctx.fillStyle = b.color;
    ctx.beginPath(); ctx.arc(p[0], p[1], rad, 0, 7); ctx.fill();
    ctx.fillStyle = "rgba(255,255,255,.85)";
    ctx.beginPath(); ctx.arc(p[0] - rad * 0.3, p[1] - rad * 0.3, rad * 0.35, 0, 7); ctx.fill();

    if (b === selected) {
      ctx.strokeStyle = "#fff";
      ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.arc(p[0], p[1], rad + 7, 0, 7); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = "#fff";
      ctx.font = "11px ui-monospace, monospace";
      ctx.fillText(b.name, p[0] + rad + 11, p[1] - 4);
    }
  }
}

function rayPath(trail) {
  ctx.beginPath();
  const n = trail.length / 2;
  for (let k = 0; k < n; k++) {
    const [vx, vy] = viewXY(trail[k * 2], trail[k * 2 + 1]);
    const p = V.sheet ? project(vx, vy, sheetZ(gridPhiAt(vx, vy))) : project(vx, vy, 0);
    k ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]);
  }
  ctx.stroke();
}

function drawPhotons() {
  ctx.lineWidth = 1.3;
  ctx.strokeStyle = "#fff7c0";
  for (let i = ghostRays.length - 1; i >= 0; i--) {
    const g = ghostRays[i];
    ctx.globalAlpha = g.alpha;
    rayPath(g.trail);
    if (!paused) g.alpha *= 0.992;
    if (g.alpha < 0.04) ghostRays.splice(i, 1);
  }
  ctx.globalAlpha = 1;
  for (const ph of photons) {
    if (ph.trail.length > 2) {
      ctx.strokeStyle = "#fff7c0";
      ctx.globalAlpha = 0.55;
      rayPath(ph.trail);
      ctx.globalAlpha = 1;
    }
    const [vx, vy] = viewXY(ph.x, ph.y);
    const p = V.sheet ? project(vx, vy, sheetZ(gridPhiAt(vx, vy))) : project(vx, vy, 0);
    ctx.fillStyle = "#fffde8";
    ctx.beginPath(); ctx.arc(p[0], p[1], 1.6, 0, 7); ctx.fill();
  }
}

function drawLaunchPreview() {
  if (!drag.active || drag.mode !== "launch") return;
  const a = project(drag.wx, drag.wy, 0), b = [drag.sx, drag.sy];
  ctx.strokeStyle = "#5ce8c5";
  ctx.lineWidth = 2;
  ctx.setLineDash([6, 5]);
  ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#5ce8c5";
  ctx.beginPath(); ctx.arc(a[0], a[1], 5, 0, 7); ctx.fill();
  const [wx2, wy2] = unproject(b[0], b[1]);
  const [lvx, lvy] = launchVelocity(drag.wx, drag.wy, wx2, wy2);
  ctx.font = "11px ui-monospace, monospace";
  ctx.fillText(`v = ${(Math.hypot(lvx, lvy) * U.KMS_PER_AUDAY).toFixed(1)} km/s`, b[0] + 10, b[1]);
}

function drawStars() {
  if (!drawStars.f || drawStars.w !== W || drawStars.h !== H) {
    drawStars.f = [];
    drawStars.w = W; drawStars.h = H;
    const rnd = mulberry32(7);
    for (let i = 0; i < 140; i++) {
      drawStars.f.push([rnd() * W, rnd() * H, rnd() * 1.3 + 0.3, rnd()]);
    }
  }
  for (const [x, y, r, tw] of drawStars.f) {
    ctx.globalAlpha = 0.22 + 0.4 * Math.abs(Math.sin(performance.now() * 5e-4 + tw * 9));
    ctx.fillStyle = "#cdd6ff";
    ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.fill();
  }
  ctx.globalAlpha = 1;
}

// ---------------- graphs ----------------

const gE = $("gEnergy").getContext("2d");
const gQ = $("gCons").getContext("2d");
const gC = $("gCurv").getContext("2d");
const gK = $("gKepler").getContext("2d");

function drawEnergyGraph() {
  const w = 300, h = 100;
  gE.clearRect(0, 0, w, h);
  if (energyHist.length < 2) return;
  let lo = Infinity, hi = -Infinity;
  for (const s of energyHist) {
    lo = Math.min(lo, s.ke, s.pe, s.e);
    hi = Math.max(hi, s.ke, s.pe, s.e);
  }
  const pad = (hi - lo) * 0.1 + 1e-300; lo -= pad; hi += pad;
  const X = i => i / (HIST_MAX - 1) * w;
  const Y = v => h - (v - lo) / (hi - lo) * h;
  for (const [key, col, lw] of [["ke", "#5ce8c5", 1.2], ["pe", "#ff6b8a", 1.2], ["e", "#ffd166", 2]]) {
    gE.strokeStyle = col; gE.lineWidth = lw;
    gE.beginPath();
    energyHist.forEach((s, i) => i ? gE.lineTo(X(i), Y(s[key])) : gE.moveTo(X(i), Y(s[key])));
    gE.stroke();
  }
}

function drawConsGraph() {
  const w = 300, h = 100;
  gQ.clearRect(0, 0, w, h);
  if (consHist.length < 2) return;
  // log10 |relative drift|, range 10⁻¹⁶ .. 10⁰
  const X = i => i / (HIST_MAX - 1) * w;
  const Y = v => {
    const lg = Math.log10(Math.abs(v) + 1e-16);
    return Math.max(2, Math.min(h - 2, h - (lg + 16) / 16 * h));
  };
  for (const [key, col] of [["dE", "#ffd166"], ["dL", "#7c9bff"], ["dP", "#5ce8c5"]]) {
    gQ.strokeStyle = col; gQ.lineWidth = 1.3;
    gQ.beginPath();
    consHist.forEach((s, i) => i ? gQ.lineTo(X(i), Y(s[key])) : gQ.moveTo(X(i), Y(s[key])));
    gQ.stroke();
  }
  gQ.fillStyle = "rgba(255,255,255,.4)";
  gQ.font = "9px ui-monospace, monospace";
  for (const lg of [-4, -8, -12]) {
    const y = 100 - (lg + 16) / 16 * 100;
    gQ.fillRect(0, y, w, 0.5);
    gQ.fillText("1e" + lg, 3, y - 2);
  }
}

function drawCurvGraph() {
  const w = 300, h = 100;
  gC.clearRect(0, 0, w, h);
  if (!bodies.length) return;
  const span = viewSpan();
  const sel = selected ? viewXY(selected.x, selected.y) : null;
  const yLine = sel ? sel[1] : V.cy;
  const NS = 120;
  const phis = [], defs = [];
  let lo = 0, dMax = 1e-300;
  for (let i = 0; i <= NS; i++) {
    const x = V.cx - span + 2 * span * i / NS;
    const phi = gridPhiAt(x, yLine);
    const [ix, iy] = frameInv(x, yLine);
    const def = 1 - clockRate(potentialAt(ix, iy), 0);    // 1 − dτ/dt
    phis.push(phi); defs.push(def);
    lo = Math.min(lo, phi);
    dMax = Math.max(dMax, def);
  }
  lo = lo || -1;
  const X = i => i / NS * w;
  gC.strokeStyle = "#7c9bff"; gC.lineWidth = 1.6;
  gC.beginPath();
  phis.forEach((phi, i) => {
    const y = 6 + (phi / lo) * (h - 16);
    i ? gC.lineTo(X(i), y) : gC.moveTo(X(i), y);
  });
  gC.stroke();
  gC.strokeStyle = "#ff9d5c"; gC.lineWidth = 1.3;
  gC.beginPath();
  defs.forEach((d, i) => {
    const y = h - 4 - (d / dMax) * (h - 12);
    i ? gC.lineTo(X(i), y) : gC.moveTo(X(i), y);
  });
  gC.stroke();
  gC.fillStyle = "rgba(255,255,255,.45)";
  gC.font = "9px ui-monospace, monospace";
  gC.fillText(`max(1−dτ/dt) = ${dMax.toExponential(1)}`, 4, 10);
  for (const b of bodies) {
    const [vx] = viewXY(b.x, b.y);
    const fx = (vx - (V.cx - span)) / (2 * span);
    if (fx >= 0 && fx <= 1) { gC.fillStyle = b.color; gC.fillRect(fx * w - 1, h - 5, 2, 5); }
  }
}

function drawKeplerGraph() {
  const w = 300, h = 100;
  gK.clearRect(0, 0, w, h);
  const c = heaviest();
  if (!c) return;
  const pts = keplerPoints.filter(p => p.a > 0 && p.T > 0);
  let aLo = 0.1, aHi = 40;
  if (pts.length) {
    aLo = Math.min(...pts.map(p => p.a)) / 2;
    aHi = Math.max(...pts.map(p => p.a)) * 2;
  }
  const mu = Gv() * c.m;
  const Tof = a => 2 * Math.PI * Math.sqrt(a * a * a / mu);
  const tLo = Tof(aLo) / 2, tHi = Tof(aHi) * 2;
  const X = a => (Math.log10(a) - Math.log10(aLo)) / (Math.log10(aHi) - Math.log10(aLo)) * w;
  const Y = T => h - (Math.log10(T) - Math.log10(tLo)) / (Math.log10(tHi) - Math.log10(tLo)) * h;
  // theory: T = 2π√(a³/GM) — a straight slope-3/2 line in log-log
  gK.strokeStyle = "rgba(255,255,255,.35)";
  gK.setLineDash([4, 4]);
  gK.beginPath();
  gK.moveTo(X(aLo), Y(Tof(aLo)));
  gK.lineTo(X(aHi), Y(Tof(aHi)));
  gK.stroke();
  gK.setLineDash([]);
  for (const p of pts) {
    gK.fillStyle = p.color;
    gK.beginPath(); gK.arc(X(p.a), Y(p.T), 3, 0, 7); gK.fill();
  }
  gK.fillStyle = "rgba(255,255,255,.45)";
  gK.font = "9px ui-monospace, monospace";
  gK.fillText(pts.length ? `${pts.length} measured orbit(s) vs T=2π√(a³/GM)` : "measures periods as orbits complete…", 4, 10);
}

// ---------------- UI ----------------

function $(id) { return document.getElementById(id); }
const fmt = v => !isFinite(v) ? "—"
  : Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(2)
  : Math.abs(v) >= 1000 ? v.toFixed(0)
  : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(3);

function syncUI() {
  $("sG").value = Math.log10(P.Gmult);
  $("sC").value = Math.log10(P.cMult);
  $("sN").value = P.n;
  $("sEps").value = P.eps;
  $("sDt").value = Math.log10(P.dt);
  $("sSub").value = P.sub;
  $("sWarp").value = V.warp;
  $("sM").value = Math.log10(V.newMass);
  $("integrator").value = P.integrator;
  $("cAdaptive").checked = P.adaptive;
  $("cPN1").checked = P.pn1;
  $("cPW").checked = P.pw;
  $("cMerge").checked = P.merge;
  $("cRot").checked = V.rotFrame;
  $("cGrid").checked = V.grid;
  $("cTrail").checked = V.trails;
  $("cGlow").checked = V.glow;
  $("cFollow").checked = V.follow;
  updateOutputs();
}

function updateOutputs() {
  $("oG").textContent = fmt(P.Gmult) + "×G";
  $("oC").textContent = fmt(P.cMult) + "×c";
  $("oN").textContent = P.n.toFixed(2);
  $("oEps").textContent = P.eps.toFixed(4) + " AU";
  $("oDt").textContent = fmt(P.dt) + " d";
  $("oSub").textContent = P.sub;
  $("oWarp").textContent = V.warp.toFixed(2);
  $("oM").textContent = fmt(V.newMass) + " M☉";

  const nTxt = P.n.toFixed(2).replace(/\.?0+$/, "");
  let force = `<b>F</b> = G·m₁·m₂ / r<sup>${nTxt}</sup> &nbsp; G = ${fmt(P.Gmult)}×6.674·10⁻¹¹ m³kg⁻¹s⁻²`;
  if (P.pn1) force += ` &nbsp;<span class="tag">+ 1PN (EIH)</span>`;
  if (P.pw) force += ` &nbsp;<span class="tag">heaviest: Φ = −GM/(r−rₛ)</span>`;
  $("eqForce").innerHTML = force;
  $("eqPot").innerHTML = Math.abs(P.n - 1) < 1e-9
    ? `<b>Φ</b>(r) = G·m·ln r`
    : `<b>Φ</b>(r) = −G·m / ${P.n === 2 ? "r" : `(${(P.n - 1).toFixed(2)}·r<sup>${(P.n - 1).toFixed(2)}</sup>)`}`;
  $("eqDil").innerHTML =
    `<b>dτ/dt</b> = √(1 + 2Φ/c² − v²/c²) &nbsp; c = ${fmt(P.cMult * 299792.458)} km/s`;
}

function bind(id, fn) {
  $(id).addEventListener("input", e => { fn(parseFloat(e.target.value)); updateOutputs(); });
}
bind("sG",   v => { P.Gmult = Math.pow(10, v); computeAccels(); resetBaselines(); });
bind("sC",   v => { P.cMult = Math.pow(10, v); computeAccels(); });
bind("sN",   v => { P.n = v; computeAccels(); resetBaselines(); });
bind("sEps", v => { P.eps = v; computeAccels(); resetBaselines(); });
bind("sDt",  v => { P.dt = Math.pow(10, v); });
bind("sSub", v => { P.sub = Math.round(v); });
bind("sWarp",v => { V.warp = v; });
bind("sM",   v => { V.newMass = Math.pow(10, v); });

$("integrator").onchange = e => { P.integrator = e.target.value; resetBaselines(); };
$("cAdaptive").onchange = e => P.adaptive = e.target.checked;
$("cPN1").onchange = e => { P.pn1 = e.target.checked; computeAccels(); resetBaselines(); updateOutputs(); };
$("cPW").onchange  = e => { P.pw = e.target.checked; computeAccels(); resetBaselines(); updateOutputs(); };
$("cMerge").onchange = e => P.merge = e.target.checked;
$("cRot").onchange = e => V.rotFrame = e.target.checked;
$("cGrid").onchange  = e => V.grid = e.target.checked;
$("cTrail").onchange = e => { V.trails = e.target.checked; if (!V.trails) bodies.forEach(b => b.trail.length = 0); };
$("cGlow").onchange  = e => V.glow = e.target.checked;
$("cFollow").onchange = e => V.follow = e.target.checked;

$("preset").onchange = e => loadPreset(e.target.value);

const btnPlay = $("btnPlay");
function setPaused(p) { paused = p; btnPlay.textContent = paused ? "▶ Play" : "⏸ Pause"; }
btnPlay.onclick = () => setPaused(!paused);
$("btnStep").onclick = () => { setPaused(true); stepFrame(); };
$("btnReset").onclick = () => loadPreset($("preset").value);
$("btnView").onclick = e => {
  V.sheet = !V.sheet;
  e.target.classList.toggle("on", V.sheet);
  e.target.textContent = V.sheet ? "🗺 Top view" : "🕳 Sheet view";
};
$("btnLaunch").onclick = e => {
  V.launch = !V.launch;
  e.target.classList.toggle("on", V.launch);
  cv.style.cursor = V.launch ? "cell" : "crosshair";
};
$("btnLight").onclick = () => {
  const span = viewSpan();
  const [ox, oy] = frameInv(V.cx - span * 1.05, V.cy);
  const f = V.rotFrame && curFrame ? curFrame.frame.theta : 0;
  emitPhotonFan(ox, oy, Math.cos(f), Math.sin(f), 0, 24, span * 1.2);
};

// ---- body editor ----

function fillEditor() {
  const b = selected;
  $("editor").style.display = b ? "" : "none";
  if (!b) return;
  $("eName").textContent = b.name;
  if (document.activeElement?.dataset?.bodyfield) return;   // don't clobber typing
  $("eM").value = b.m.toExponential(4);
  $("eX").value = b.x.toPrecision(6);
  $("eY").value = b.y.toPrecision(6);
  $("eVx").value = b.vx.toPrecision(6);
  $("eVy").value = b.vy.toPrecision(6);
}
$("btnApplyBody").onclick = () => {
  const b = selected;
  if (!b) return;
  const g = id => parseFloat($(id).value);
  const vals = [g("eM"), g("eX"), g("eY"), g("eVx"), g("eVy")];
  if (vals.some(v => !isFinite(v))) return;
  [b.m, b.x, b.y, b.vx, b.vy] = vals;
  b.radAU = bodyRadiusAU(b.m);
  b.trail.length = 0; b.peri.length = 0; b._angPrev = null; b._pr1 = b._pr2 = null;
  computeAccels(); resetBaselines();
};
$("btnDelBody").onclick = () => {
  if (!selected) return;
  bodies = bodies.filter(b => b !== selected);
  selected = null;
  computeAccels(); resetBaselines();
};

// ---- export / import ----

function download(name, text, type) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name;
  a.click();
  URL.revokeObjectURL(a.href);
}
$("btnCSV").onclick = () => {
  const ids = bodies.map(b => [b.id, b.name.replace(/[^\w]+/g, "_")]);
  let csv = "t_days," + ids.map(([, n]) => `${n}.x_AU,${n}.y_AU,${n}.vx_AUday,${n}.vy_AUday`).join(",") + "\n";
  for (const row of recording) {
    const by = {};
    for (const s of row.s) by[s[0]] = s;
    csv += row.t.toPrecision(9) + "," + ids.map(([id]) =>
      by[id] ? by[id].slice(1).map(v => v.toPrecision(9)).join(",") : ",,,").join(",") + "\n";
  }
  download("gravity-lab-trajectories.csv", csv, "text/csv");
};
$("btnSave").onclick = () => {
  download("gravity-lab-state.json", JSON.stringify({
    version: 2, rngSeed, time,
    params: { ...P },
    view: { scale: V.scale, cx: V.cx, cy: V.cy, sheet: V.sheet, rotFrame: V.rotFrame },
    bodies: bodies.map(b => ({ name: b.name, m: b.m, x: b.x, y: b.y, vx: b.vx, vy: b.vy, color: b.color })),
  }, null, 1), "application/json");
};
$("fileLoad").onchange = e => {
  const f = e.target.files[0];
  if (!f) return;
  f.text().then(txt => {
    const d = JSON.parse(txt);
    bodies = []; photons = []; selected = null; colorIdx = 0; bodyIdSeq = 1;
    keplerPoints.length = 0; ghostRays.length = 0;
    recording = []; recStride = 1; recCount = 0;
    time = d.time || 0;
    Object.assign(P, d.params || {});
    Object.assign(V, d.view || {});
    for (const s of d.bodies || []) {
      bodies.push(makeBody(s.m, s.x, s.y, s.vx, s.vy, s.color, s.name));
    }
    computeAccels(); resetBaselines(); syncUI();
    $("presetInfo").textContent = `Loaded ${f.name} (${bodies.length} bodies)`;
  }).catch(err => { $("presetInfo").textContent = "Load failed: " + err.message; });
  e.target.value = "";
};
$("btnLoad").onclick = () => $("fileLoad").click();

window.addEventListener("keydown", e => {
  if (e.code === "Space" && !["INPUT", "SELECT"].includes(e.target.tagName)) {
    e.preventDefault(); setPaused(!paused);
  }
});

// ---------------- mouse ----------------

function launchVelocity(wx0, wy0, wx1, wy1) {   // view-space drag → inertial velocity
  const span = Math.max(viewSpan(), 1e-9);
  let Mtot = 0;
  for (const b of bodies) Mtot += b.m;
  const vChar = Mtot > 0 ? Math.sqrt(Gv() * Mtot / span) : Math.sqrt(Gv() / span);
  const k = vChar * 3 / span;
  let vx = (wx1 - wx0) * k, vy = (wy1 - wy0) * k;
  if (V.rotFrame && curFrame) {
    // rotate into inertial axes and add the frame's co-rotation velocity
    const f = curFrame.frame;
    const ct = Math.cos(f.theta), st = Math.sin(f.theta);
    [vx, vy] = [vx * ct - vy * st, vx * st + vy * ct];
    const [ix, iy] = frameInv(wx0, wy0);
    vx += -f.omega * (iy - f.by);
    vy += f.omega * (ix - f.bx);
  }
  return [vx, vy];
}

const drag = { active: false, mode: "", sx: 0, sy: 0, wx: 0, wy: 0 };

cv.addEventListener("mousedown", e => {
  const [wx, wy] = unproject(e.offsetX, e.offsetY);
  let hit = null;
  for (const b of bodies) {
    const [vx, vy] = viewXY(b.x, b.y);
    const p = V.sheet ? project(vx, vy, sheetZ(gridPhiAt(vx, vy))) : project(vx, vy, 0);
    const rad = bodyDrawRadius(b) + 6;
    if ((p[0] - e.offsetX) ** 2 + (p[1] - e.offsetY) ** 2 < rad * rad) hit = b;
  }
  if (hit && !V.launch) { selected = hit; return; }

  drag.active = true;
  drag.sx = e.offsetX; drag.sy = e.offsetY;
  drag.wx = wx; drag.wy = wy;
  drag.mode = V.launch ? "launch" : (V.sheet ? "orbit" : "pan");
  if (!hit && !V.launch) selected = null;
});

cv.addEventListener("mousemove", e => {
  if (!drag.active) return;
  const dx = e.offsetX - drag.sx, dy = e.offsetY - drag.sy;
  if (drag.mode === "pan") {
    V.cx -= dx / V.scale; V.cy -= dy / V.scale;
    drag.sx = e.offsetX; drag.sy = e.offsetY;
  } else if (drag.mode === "orbit") {
    V.yaw += dx * 0.005;
    V.pitch = Math.max(0.25, Math.min(1.45, V.pitch + dy * 0.005));
    drag.sx = e.offsetX; drag.sy = e.offsetY;
  } else if (drag.mode === "launch") {
    drag.sx = e.offsetX; drag.sy = e.offsetY;
  }
});

window.addEventListener("mouseup", () => {
  if (drag.active && drag.mode === "launch") {
    const [wx2, wy2] = unproject(drag.sx, drag.sy);
    const [vx, vy] = launchVelocity(drag.wx, drag.wy, wx2, wy2);
    const [ix, iy] = frameInv(drag.wx, drag.wy);
    const b = makeBody(V.newMass, ix, iy, vx, vy);
    bodies.push(b);
    selected = b;
    computeAccels(); resetBaselines();
  }
  drag.active = false;
});

cv.addEventListener("wheel", e => {
  e.preventDefault();
  const k = Math.exp(-e.deltaY * 0.0012);
  if (!V.sheet) {
    const [wx, wy] = unproject(e.offsetX, e.offsetY);
    V.scale *= k;
    const [wx2, wy2] = unproject(e.offsetX, e.offsetY);
    V.cx += wx - wx2; V.cy += wy - wy2;
  } else {
    V.scale *= k;
  }
  V.scale = Math.max(0.001, Math.min(50000, V.scale));
}, { passive: false });

// ---------------- diagnostics text ----------------

function timeText(d) {
  return d < 2 * U.YEAR_D ? `${d.toFixed(1)} d` : `${d.toFixed(0)} d (${(d / U.YEAR_D).toFixed(2)} yr)`;
}

function updateReadout() {
  const en = energyHist.length ? energyHist[energyHist.length - 1] : { ke: 0, pe: 0, e: 0 };
  const cs = consHist.length ? consHist[consHist.length - 1] : { dE: 0, dL: 0, dP: 0 };
  let txt =
    `t        = ${timeText(time)}\n` +
    `bodies   = ${bodies.length}${photons.length ? `  (+${photons.length} photons)` : ""}\n` +
    `E        = ${en.e.toExponential(3)} M☉AU²/d²\n` +
    `ΔE/E₀    = ${cs.dE.toExponential(2)}${P.pn1 ? " (Newtonian E; 1PN adds terms)" : ""}\n` +
    `ΔL/L₀    = ${cs.dL.toExponential(2)}\n` +
    `|p|/p̂    = ${cs.dP.toExponential(2)}\n` +
    `substeps = ${P._steps || 1}/step (${P.integrator}${P.adaptive ? ", adaptive" : ""})`;
  $("readout").textContent = txt;

  // selected-body science
  const b = selected;
  let s = "";
  if (b) {
    const c = heaviest() === b ? null : heaviest();
    const phi = potentialAt(b.x, b.y, b);
    const v2 = b.vx * b.vx + b.vy * b.vy;
    const rate = clockRate(phi, v2);
    s += `${b.name}\n`;
    s += `m  = ${b.m.toExponential(3)} M☉ = ${(b.m * U.MSUN_KG).toExponential(3)} kg\n`;
    s += `|v| = ${(Math.sqrt(v2) * U.KMS_PER_AUDAY).toFixed(2)} km/s\n`;
    s += `1−dτ/dt = ${(1 - rate).toExponential(3)}\n`;
    if (c) {
      const el = orbitalElements(b, c);
      if (el && el.bound) {
        s += `a  = ${fmt(el.a)} AU   e = ${el.e.toFixed(4)}\n`;
        s += `T  = ${timeText(el.T)} (Kepler)\n`;
        if (b.Tmeas) s += `T  = ${timeText(b.Tmeas)} (measured)\n`;
        const meas = measuredPrecession(b);
        const pred = predictedPrecession(b, c);
        if (meas !== null) s += `Δϖ = ${(meas * U.ARCSEC).toFixed(2)}″/orbit (measured, ${b.peri.length} passages)\n`;
        if (pred !== null) s += `Δϖ = ${(pred * U.ARCSEC).toFixed(2)}″/orbit (1PN prediction)\n`;
      } else {
        s += `unbound w.r.t. ${c.name} (hyperbolic)\n`;
      }
    }
  }
  $("selInfo").textContent = s;
  fillEditor();
}

// ---------------- main loop ----------------

let frameCount = 0;

function frame() {
  curFrame = (V.rotFrame && bodies.length >= 2) ? lagrangeData() : null;
  if (!paused && bodies.length) stepFrame();
  if (V.follow && selected) {
    const [vx, vy] = viewXY(selected.x, selected.y);
    V.cx = vx; V.cy = vy;
  }

  ctx.fillStyle = "#05060e";
  ctx.fillRect(0, 0, W, H);
  drawStars();
  drawGrid();
  drawBHmarkers();
  drawLagrangePoints();
  drawBodies();
  drawPhotons();
  drawLaunchPreview();

  if (frameCount++ % 3 === 0) {
    updateReadout();
    drawEnergyGraph();
    drawConsGraph();
    drawCurvGraph();
    drawKeplerGraph();
  }
  requestAnimationFrame(frame);
}

loadPreset("solar");
updateOutputs();
requestAnimationFrame(frame);
