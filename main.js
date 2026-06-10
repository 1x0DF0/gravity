/* ============================================================
   Gravity Lab — N-body spacetime curvature simulator
   Physics: velocity-Verlet integration of
       F = G·m1·m2 / r_s^n   (Plummer-softened, r_s = √(r²+ε²))
   Visualization is driven directly by the Newtonian potential
       Φ(p) = Σ −G·m / ((n−1)·r_s^(n−1))      (n ≠ 1)
   Weak-field gravitational time dilation:
       dτ/dt = √(1 + 2Φ/c² − v²/c²)
   ============================================================ */

"use strict";

// ---------------- state ----------------

const P = {                 // physics parameters (user-editable equations)
  G: 1,
  n: 2,                     // force exponent: F ∝ 1/rⁿ
  eps: 2,                   // softening length ε
  dt: 0.05,
  sub: 8,                   // integration substeps per frame
  c: 30,                    // speed of light (sim units) for dilation display
  merge: false,
};

const V = {                 // view / visual parameters
  sheet: false,             // false = top-down 2D, true = 3D embedding sheet
  warp: 1,
  grid: true,
  trails: true,
  glow: true,
  follow: false,
  yaw: 0.6, pitch: 1.05,    // 3D sheet orientation
  scale: 1.6,               // pixels per world unit
  cx: 0, cy: 0,             // camera center (world coords)
  launch: false,
  newMass: 10,
};

let bodies = [];
let time = 0;
let paused = false;
let selected = null;
let E0 = null;              // initial total energy (for drift readout)
const energyHist = [];      // {ke, pe, e}
const HIST_MAX = 400;
const TRAIL_MAX = 350;

const PALETTE = ["#ffd166","#7c9bff","#ff6b8a","#5ce8c5","#c792ea","#ff9d5c","#7fd4ff","#aef07a"];
let colorIdx = 0;

function makeBody(m, x, y, vx, vy, color) {
  return {
    m, x, y, vx, vy,
    ax: 0, ay: 0,
    r: Math.max(2.2, 2.2 * Math.cbrt(m)),
    color: color || PALETTE[colorIdx++ % PALETTE.length],
    trail: [],
  };
}

// ---------------- physics (exact to the equations shown) ----------------

function computeAccels() {
  const N = bodies.length;
  for (const b of bodies) { b.ax = 0; b.ay = 0; }
  const e2 = P.eps * P.eps;
  for (let i = 0; i < N; i++) {
    const bi = bodies[i];
    for (let j = i + 1; j < N; j++) {
      const bj = bodies[j];
      const dx = bj.x - bi.x, dy = bj.y - bi.y;
      const rs2 = dx * dx + dy * dy + e2;
      // a_i = G·m_j·(r_vec) / r_s^(n+1)  — reduces to Plummer softening at n=2
      let inv;
      if (P.n === 2) { const rs = Math.sqrt(rs2); inv = 1 / (rs2 * rs); }
      else inv = Math.pow(rs2, -(P.n + 1) / 2);
      const f = P.G * inv;
      bi.ax += f * bj.m * dx;  bi.ay += f * bj.m * dy;
      bj.ax -= f * bi.m * dx;  bj.ay -= f * bi.m * dy;
    }
  }
}

// pair potential consistent with the softened force (so energy is conserved)
function pairPotential(mi, mj, rs) {
  if (Math.abs(P.n - 1) < 1e-9) return P.G * mi * mj * Math.log(rs);
  return -P.G * mi * mj / ((P.n - 1) * Math.pow(rs, P.n - 1));
}

function energies() {
  let ke = 0, pe = 0;
  const e2 = P.eps * P.eps;
  for (let i = 0; i < bodies.length; i++) {
    const bi = bodies[i];
    ke += 0.5 * bi.m * (bi.vx * bi.vx + bi.vy * bi.vy);
    for (let j = i + 1; j < bodies.length; j++) {
      const bj = bodies[j];
      const dx = bj.x - bi.x, dy = bj.y - bi.y;
      pe += pairPotential(bi.m, bj.m, Math.sqrt(dx * dx + dy * dy + e2));
    }
  }
  return { ke, pe, e: ke + pe };
}

// gravitational potential per unit mass at a point (drives all curvature visuals)
// pass `skip` to exclude a body's own (infinite/self) contribution at its location
function potentialAt(x, y, skip) {
  const e2 = P.eps * P.eps;
  let phi = 0;
  for (const b of bodies) {
    if (b === skip) continue;
    const dx = b.x - x, dy = b.y - y;
    const rs = Math.sqrt(dx * dx + dy * dy + e2);
    phi += (Math.abs(P.n - 1) < 1e-9)
      ? P.G * b.m * Math.log(rs)
      : -P.G * b.m / ((P.n - 1) * Math.pow(rs, P.n - 1));
  }
  return phi;
}

function clockRate(phi, v2) {           // dτ/dt, clamped at the "horizon"
  return Math.sqrt(Math.max(0, 1 + 2 * phi / (P.c * P.c) - v2 / (P.c * P.c)));
}

function verletStep(dt) {
  for (const b of bodies) {
    b.x += b.vx * dt + 0.5 * b.ax * dt * dt;
    b.y += b.vy * dt + 0.5 * b.ay * dt * dt;
    b.ox = b.ax; b.oy = b.ay;
  }
  computeAccels();
  for (const b of bodies) {
    b.vx += 0.5 * (b.ox + b.ax) * dt;
    b.vy += 0.5 * (b.oy + b.ay) * dt;
  }
  time += dt;
}

function mergeCollisions() {
  for (let i = 0; i < bodies.length; i++) {
    for (let j = i + 1; j < bodies.length; j++) {
      const a = bodies[i], b = bodies[j];
      const dx = b.x - a.x, dy = b.y - a.y;
      if (dx * dx + dy * dy < (a.r + b.r) * (a.r + b.r) * 0.36) {
        const m = a.m + b.m;                       // conserve momentum
        a.x = (a.x * a.m + b.x * b.m) / m;
        a.y = (a.y * a.m + b.y * b.m) / m;
        a.vx = (a.vx * a.m + b.vx * b.m) / m;
        a.vy = (a.vy * a.m + b.vy * b.m) / m;
        a.m = m;
        a.r = Math.max(2.2, 2.2 * Math.cbrt(m));
        if (selected === b) selected = a;
        bodies.splice(j, 1); j--;
      }
    }
  }
}

function stepFrame() {
  for (let s = 0; s < P.sub; s++) verletStep(P.dt);
  if (P.merge) mergeCollisions();
  for (const b of bodies) {
    // record Φ at the time of the sample so 3D trails ride the historical sheet
    b.trail.push(b.x, b.y, potentialAt(b.x, b.y, b));
    if (b.trail.length > TRAIL_MAX * 3) b.trail.splice(0, b.trail.length - TRAIL_MAX * 3);
  }
  const en = energies();
  if (E0 === null && bodies.length) E0 = en.e;
  energyHist.push(en);
  if (energyHist.length > HIST_MAX) energyHist.shift();
}

// ---------------- presets ----------------

function circV(M, r) { return Math.sqrt(P.G * M / r); }   // circular-orbit speed (n=2)

const PRESETS = {
  orbit() {
    P.G = 1; P.n = 2; P.eps = 2; P.dt = 0.05; V.scale = 1.6;
    const M = 1200;
    bodies.push(makeBody(M, 0, 0, 0, 0, "#ffd166"));
    for (const [r, m, col] of [[120, 4, "#7c9bff"], [210, 9, "#5ce8c5"]]) {
      bodies.push(makeBody(m, r, 0, 0, circV(M, r), col));
    }
  },
  binary() {
    P.G = 1; P.n = 2; P.eps = 2; P.dt = 0.04; V.scale = 2.2;
    const m = 500, d = 90;
    const v = Math.sqrt(P.G * m / (2 * d));   // circular about the barycenter
    bodies.push(makeBody(m, -d / 2, 0, 0, -v, "#ffd166"));
    bodies.push(makeBody(m,  d / 2, 0, 0,  v, "#ff6b8a"));
    bodies.push(makeBody(2, 0, 200, circV(2 * m, 200), 0, "#7fd4ff"));
  },
  figure8() {
    // Chenciner–Montgomery figure-eight choreography (G = 1, m = 1, exact ICs)
    P.G = 1; P.n = 2; P.eps = 0.02; P.dt = 0.004; P.sub = 20; V.scale = 240;
    const p = [[-0.97000436, 0.24308753], [0.97000436, -0.24308753], [0, 0]];
    const v3 = [-0.93240737, -0.86473146];
    const vel = [[-v3[0] / 2, -v3[1] / 2], [-v3[0] / 2, -v3[1] / 2], v3];
    for (let i = 0; i < 3; i++) {
      bodies.push(makeBody(1, p[i][0], p[i][1], vel[i][0], vel[i][1]));
    }
  },
  solar() {
    P.G = 1; P.n = 2; P.eps = 1.5; P.dt = 0.04; V.scale = 1.1;
    const M = 2500;
    bodies.push(makeBody(M, 0, 0, 0, 0, "#ffd166"));
    const orbits = [[70, 1.5], [115, 3], [170, 5], [240, 8], [330, 16]];
    orbits.forEach(([r, m], i) => {
      const a = i * 2.39996;                     // golden-angle spread
      const v = circV(M, r);
      bodies.push(makeBody(m, r * Math.cos(a), r * Math.sin(a), -v * Math.sin(a), v * Math.cos(a)));
    });
  },
  slingshot() {
    P.G = 1; P.n = 2; P.eps = 3; P.dt = 0.03; V.scale = 1.4;
    bodies.push(makeBody(4000, 0, 0, 0, 0, "#1a1a2e"));
    bodies.push(makeBody(6, -380, 90, 9, 0, "#5ce8c5"));
    bodies.push(makeBody(6, -380, -55, 8, 0, "#ff9d5c"));
    bodies.push(makeBody(40, 250, 0, 0, circV(4000, 250), "#7c9bff"));
  },
  chaos() {
    P.G = 1; P.n = 2; P.eps = 4; P.dt = 0.03; V.scale = 1.3;
    for (let i = 0; i < 12; i++) {
      const a = Math.random() * Math.PI * 2, r = 60 + Math.random() * 220;
      const m = 5 + Math.random() * 120;
      const v = circV(900, r) * (0.6 + Math.random() * 0.5);
      bodies.push(makeBody(m, r * Math.cos(a), r * Math.sin(a), -v * Math.sin(a), v * Math.cos(a)));
    }
  },
  empty() { P.G = 1; P.n = 2; V.scale = 1.6; },
};

function loadPreset(name) {
  bodies = []; colorIdx = 0; time = 0; E0 = null;
  energyHist.length = 0; selected = null;
  V.cx = 0; V.cy = 0; P.sub = 8;
  (PRESETS[name] || PRESETS.orbit)();
  computeAccels();
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

const VIEW_W = () => W - 340;   // canvas area left of the panel

// world → screen, always [sx, sy, depth].  In sheet view, z comes from the potential field.
function project(x, y, z) {
  const dx = x - V.cx, dy = y - V.cy;
  if (!V.sheet) return [VIEW_W() / 2 + dx * V.scale, H / 2 + dy * V.scale, 0];
  const cy_ = Math.cos(V.yaw), sy_ = Math.sin(V.yaw);
  const xr = dx * cy_ - dy * sy_;
  const yr = dx * sy_ + dy * cy_;
  const cp = Math.cos(V.pitch), sp = Math.sin(V.pitch);
  return [
    VIEW_W() / 2 + xr * V.scale,
    H / 2 + (yr * cp - (z || 0) * sp) * V.scale,
    yr,                                          // depth, for painter's sort
  ];
}

function unproject(sx, sy) {                     // screen → world on the z=0 plane
  const px = (sx - VIEW_W() / 2) / V.scale;
  const py = (sy - H / 2) / V.scale;
  if (!V.sheet) return [px + V.cx, py + V.cy];
  const cp = Math.cos(V.pitch);
  const xr = px, yr = py / (cp || 1e-6);
  const cy_ = Math.cos(-V.yaw), sy_ = Math.sin(-V.yaw);
  return [xr * cy_ - yr * sy_ + V.cx, xr * sy_ + yr * cy_ + V.cy];
}

// Sheet depth is the potential, self-normalized to the values currently on
// screen so the visualization stays readable at any G / mass / zoom level.
// (Visual only — the physics integrates the un-normalized equations.)
let phiLo = -1, phiHi = 0, gridSpan = 200;   // refreshed each frame by drawGrid
function depthT(phi) {
  return Math.max(0, Math.min(1, (phiHi - phi) / (phiHi - phiLo + 1e-12)));
}
function sheetZ(phi) {
  // depth capped in screen pixels so deep wells never fall off the canvas
  const maxDepth = Math.min(V.warp * gridSpan * 0.45,
                            H * 0.3 / (V.scale * Math.max(0.3, Math.sin(V.pitch))));
  return -maxDepth * Math.pow(depthT(phi), 0.75);
}

function drawRadius(b) {
  return b.r * Math.min(2.0, Math.max(0.5, V.scale * 0.7));
}

// ---------------- rendering ----------------

function heatColor(t) {  // 0 = flat space, 1 = deep well
  t = Math.max(0, Math.min(1, t));
  const r = Math.round(40 + 215 * Math.pow(t, 1.4));
  const g = Math.round(60 + 90 * t * (1 - t) * 4 * 0.6);
  const b = Math.round(140 + 100 * (1 - t));
  return `rgb(${r},${g},${b})`;
}

function drawGrid() {
  const NG = V.sheet ? 36 : 30;
  const span = Math.max(VIEW_W(), H) / V.scale * (V.sheet ? 0.85 : 0.6);
  const step = (2 * span) / NG;
  const x0 = V.cx - span, y0 = V.cy - span;
  gridSpan = span;

  // pass 1: sample the potential on the lattice, find its on-screen range
  const phis = [];
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i <= NG; i++) {
    for (let j = 0; j <= NG; j++) {
      const phi = potentialAt(x0 + i * step, y0 + j * step);
      phis.push(phi);
      if (phi < lo) lo = phi;
      if (phi > hi) hi = phi;
    }
  }
  phiLo = lo; phiHi = hi;
  if (!V.grid) return;   // normalization above is still needed by the sheet/bodies

  // pass 2: place nodes — sheet view drops them by Φ, 2D view pulls them toward the masses
  const pts = [];   // [sx, sy, depth, t] per node
  for (let i = 0; i <= NG; i++) {
    for (let j = 0; j <= NG; j++) {
      const k = i * (NG + 1) + j;
      const wx = x0 + i * step, wy = y0 + j * step;
      const t = depthT(phis[k]);
      let p;
      if (V.sheet) {
        p = project(wx, wy, sheetZ(phis[k]));
      } else {
        let ox = 0, oy = 0;   // direction of the pull, weighted per body
        for (const b of bodies) {
          const dx = b.x - wx, dy = b.y - wy;
          const d = Math.sqrt(dx * dx + dy * dy) + 1e-9;
          const w = b.m / Math.pow(d + P.eps, P.n - 1);
          ox += dx / d * w; oy += dy / d * w;
        }
        const len = Math.sqrt(ox * ox + oy * oy) + 1e-12;
        const pull = V.warp * step * 1.6 * Math.pow(t, 0.9);
        p = project(wx + ox / len * pull, wy + oy / len * pull, 0);
      }
      p[3] = Math.sqrt(t);   // widened response so the well is visible far out
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

function drawBodies() {
  const order = [...bodies];
  if (V.sheet) order.sort((a, b) => project(a.x, a.y, 0)[2] - project(b.x, b.y, 0)[2]);

  for (const b of order) {
    // trail
    if (V.trails && b.trail.length > 6) {
      ctx.lineWidth = 1.5;
      const n = b.trail.length / 3;
      ctx.strokeStyle = b.color;
      ctx.beginPath();
      for (let k = 0; k < n; k++) {
        const tx = b.trail[k * 3], ty = b.trail[k * 3 + 1];
        const pp = V.sheet ? project(tx, ty, sheetZ(b.trail[k * 3 + 2])) : project(tx, ty, 0);
        k ? ctx.lineTo(pp[0], pp[1]) : ctx.moveTo(pp[0], pp[1]);
      }
      ctx.globalAlpha = 0.35;
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    const phi = potentialAt(b.x, b.y, b);
    const p = V.sheet ? project(b.x, b.y, sheetZ(phi)) : project(b.x, b.y, 0);
    const rad = drawRadius(b);

    if (V.glow) {
      const g = ctx.createRadialGradient(p[0], p[1], 0, p[0], p[1], rad * 5);
      g.addColorStop(0, b.color + "");
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
      const v2 = b.vx * b.vx + b.vy * b.vy;
      const rate = clockRate(phi, v2);
      ctx.fillStyle = "#fff";
      ctx.font = "11px ui-monospace, monospace";
      ctx.fillText(`m=${fmt(b.m)}  |v|=${fmt(Math.sqrt(v2))}  dτ/dt=${rate.toFixed(4)}`, p[0] + rad + 12, p[1] - 4);
    }
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
  const [vx, vy] = launchVelocity(drag.wx, drag.wy, wx2, wy2);
  ctx.font = "11px ui-monospace, monospace";
  ctx.fillText(`v = (${fmt(vx)}, ${fmt(vy)})`, b[0] + 10, b[1]);
}

function drawStars() {
  // static starfield, regenerated on resize
  if (!drawStars.f || drawStars.w !== W || drawStars.h !== H) {
    drawStars.f = [];
    drawStars.w = W; drawStars.h = H;
    for (let i = 0; i < 140; i++) {
      drawStars.f.push([Math.random() * W, Math.random() * H, Math.random() * 1.3 + 0.3, Math.random()]);
    }
  }
  for (const [x, y, r, tw] of drawStars.f) {
    ctx.globalAlpha = 0.25 + 0.45 * Math.abs(Math.sin(time * 0.7 + tw * 9));
    ctx.fillStyle = "#cdd6ff";
    ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.fill();
  }
  ctx.globalAlpha = 1;
}

// ---------------- graphs ----------------

const gE = document.getElementById("gEnergy").getContext("2d");
const gC = document.getElementById("gCurv").getContext("2d");

function drawEnergyGraph() {
  const w = 300, h = 110;
  gE.clearRect(0, 0, w, h);
  if (energyHist.length < 2) return;
  let lo = Infinity, hi = -Infinity;
  for (const s of energyHist) {
    lo = Math.min(lo, s.ke, s.pe, s.e);
    hi = Math.max(hi, s.ke, s.pe, s.e);
  }
  const pad = (hi - lo) * 0.1 + 1e-9; lo -= pad; hi += pad;
  const X = i => i / (HIST_MAX - 1) * w;
  const Y = v => h - (v - lo) / (hi - lo) * h;
  const series = [["ke", "#5ce8c5"], ["pe", "#ff6b8a"], ["e", "#ffd166"]];
  for (const [key, col] of series) {
    gE.strokeStyle = col;
    gE.lineWidth = key === "e" ? 2 : 1.2;
    gE.beginPath();
    energyHist.forEach((s, i) => i ? gE.lineTo(X(i), Y(s[key])) : gE.moveTo(X(i), Y(s[key])));
    gE.stroke();
  }
  if (Y(0) > 0 && Y(0) < h) {
    gE.strokeStyle = "rgba(255,255,255,.15)";
    gE.beginPath(); gE.moveTo(0, Y(0)); gE.lineTo(w, Y(0)); gE.stroke();
  }
}

function drawCurvGraph() {
  const w = 300, h = 110;
  gC.clearRect(0, 0, w, h);
  if (!bodies.length) return;
  const span = Math.max(VIEW_W(), 100) / V.scale * 0.55;
  const yLine = selected ? selected.y : V.cy;
  const NS = 120;
  const phis = [], rates = [];
  let lo = 0;
  for (let i = 0; i <= NS; i++) {
    const x = V.cx - span + (2 * span) * i / NS;
    const phi = potentialAt(x, yLine);
    phis.push(phi); rates.push(clockRate(phi, 0));
    lo = Math.min(lo, phi);
  }
  lo = lo || -1;
  const X = i => i / NS * w;
  // potential well Φ(x)
  gC.strokeStyle = "#7c9bff"; gC.lineWidth = 1.6;
  gC.beginPath();
  phis.forEach((phi, i) => {
    const y = 8 + (phi / lo) * (h - 20);
    i ? gC.lineTo(X(i), y) : gC.moveTo(X(i), y);
  });
  gC.stroke();
  // clock rate dτ/dt ∈ [0,1]
  gC.strokeStyle = "#ff9d5c"; gC.lineWidth = 1.3;
  gC.beginPath();
  rates.forEach((r, i) => {
    const y = h - 6 - r * (h - 14);
    i ? gC.lineTo(X(i), y) : gC.moveTo(X(i), y);
  });
  gC.stroke();
  // body markers
  gC.fillStyle = "rgba(255,255,255,.5)";
  for (const b of bodies) {
    const fx = (b.x - (V.cx - span)) / (2 * span);
    if (fx >= 0 && fx <= 1) gC.fillRect(fx * w - 1, h - 5, 2, 5);
  }
}

// ---------------- UI ----------------

const $ = id => document.getElementById(id);
const fmt = v => Math.abs(v) >= 1000 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2);

function syncUI() {
  $("sG").value = Math.log10(P.G);
  $("sN").value = P.n;
  $("sEps").value = P.eps;
  $("sDt").value = Math.log10(P.dt);
  $("sSub").value = P.sub;
  $("sC").value = Math.log10(P.c);
  $("sWarp").value = V.warp;
  $("sM").value = Math.log10(V.newMass);
  updateOutputs();
}

function updateOutputs() {
  $("oG").textContent = fmt(P.G);
  $("oN").textContent = P.n.toFixed(2);
  $("oEps").textContent = P.eps.toFixed(1);
  $("oDt").textContent = P.dt.toFixed(3);
  $("oSub").textContent = P.sub;
  $("oC").textContent = fmt(P.c);
  $("oWarp").textContent = V.warp.toFixed(2);
  $("oM").textContent = fmt(V.newMass);
  const nTxt = P.n.toFixed(2).replace(/\.?0+$/, "");
  $("eqForce").innerHTML = `<b>F</b> = G·m₁·m₂ / r<sup>${nTxt}</sup> &nbsp; (G = ${fmt(P.G)}, ε = ${P.eps.toFixed(1)})`;
  $("eqPot").innerHTML = Math.abs(P.n - 1) < 1e-9
    ? `<b>Φ</b>(r) = G·m·ln r`
    : `<b>Φ</b>(r) = −G·m / ${P.n === 2 ? "r" : `(${(P.n - 1).toFixed(2)}·r<sup>${(P.n - 1).toFixed(2)}</sup>)`}`;
  $("eqDil").innerHTML = `<b>dτ/dt</b> = √(1 + 2Φ/c² − v²/c²) &nbsp; (c = ${fmt(P.c)})`;
}

function bind(id, fn) { $(id).addEventListener("input", e => { fn(parseFloat(e.target.value)); updateOutputs(); }); }
bind("sG",   v => { P.G = Math.pow(10, v); computeAccels(); });
bind("sN",   v => { P.n = v; computeAccels(); E0 = null; energyHist.length = 0; });
bind("sEps", v => { P.eps = v; computeAccels(); E0 = null; energyHist.length = 0; });
bind("sDt",  v => { P.dt = Math.pow(10, v); });
bind("sSub", v => { P.sub = Math.round(v); });
bind("sC",   v => { P.c = Math.pow(10, v); });
bind("sWarp",v => { V.warp = v; });
bind("sM",   v => { V.newMass = Math.pow(10, v); });

$("cGrid").onchange  = e => V.grid = e.target.checked;
$("cTrail").onchange = e => { V.trails = e.target.checked; if (!V.trails) bodies.forEach(b => b.trail.length = 0); };
$("cGlow").onchange  = e => V.glow = e.target.checked;
$("cMerge").onchange = e => P.merge = e.target.checked;
$("cFollow").onchange= e => V.follow = e.target.checked;

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

window.addEventListener("keydown", e => {
  if (e.code === "Space" && e.target.tagName !== "INPUT" && e.target.tagName !== "SELECT") {
    e.preventDefault(); setPaused(!paused);
  }
});

// ---------------- mouse ----------------

// Drag → velocity, scaled to the system's characteristic orbital speed so
// launching feels right whether masses are 1 or 10⁶.
function launchVelocity(wx0, wy0, wx1, wy1) {
  const span = Math.max(VIEW_W() / V.scale / 2, 1e-6);
  let Mtot = 0;
  for (const b of bodies) Mtot += b.m;
  const vChar = Mtot > 0 ? Math.sqrt(P.G * Mtot / span) : 1;
  const k = vChar * 3 / span;
  return [(wx1 - wx0) * k, (wy1 - wy0) * k];
}

const drag = { active: false, mode: "", sx: 0, sy: 0, wx: 0, wy: 0 };

cv.addEventListener("mousedown", e => {
  const [wx, wy] = unproject(e.offsetX, e.offsetY);
  // body picking (screen-space)
  let hit = null;
  for (const b of bodies) {
    const p = V.sheet ? project(b.x, b.y, sheetZ(potentialAt(b.x, b.y, b))) : project(b.x, b.y, 0);
    const rad = drawRadius(b) + 6;
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
    drag.sx = e.offsetX; drag.sy = e.offsetY;   // current endpoint for preview
  }
});

window.addEventListener("mouseup", e => {
  if (drag.active && drag.mode === "launch") {
    const [wx2, wy2] = unproject(drag.sx, drag.sy);
    const [vx, vy] = launchVelocity(drag.wx, drag.wy, wx2, wy2);
    bodies.push(makeBody(V.newMass, drag.wx, drag.wy, vx, vy));
    computeAccels();
    E0 = null; energyHist.length = 0;
  }
  drag.active = false;
});

cv.addEventListener("wheel", e => {
  e.preventDefault();
  const k = Math.exp(-e.deltaY * 0.0012);
  // zoom about the cursor (2D); about center in sheet view
  if (!V.sheet) {
    const [wx, wy] = unproject(e.offsetX, e.offsetY);
    V.scale *= k;
    const [wx2, wy2] = unproject(e.offsetX, e.offsetY);
    V.cx += wx - wx2; V.cy += wy - wy2;
  } else {
    V.scale *= k;
  }
  V.scale = Math.max(0.02, Math.min(3000, V.scale));
}, { passive: false });

// ---------------- main loop ----------------

function frame() {
  if (!paused && bodies.length) stepFrame();
  if (V.follow && selected) { V.cx = selected.x; V.cy = selected.y; }

  // background
  ctx.fillStyle = "#05060e";
  ctx.fillRect(0, 0, W, H);
  drawStars();
  drawGrid();
  drawBodies();
  drawLaunchPreview();

  // readout
  const en = energyHist.length ? energyHist[energyHist.length - 1] : { ke: 0, pe: 0, e: 0 };
  const drift = (E0 !== null && Math.abs(E0) > 1e-9) ? ((en.e - E0) / Math.abs(E0) * 100) : 0;
  $("readout").textContent =
    `t      = ${time.toFixed(1)}\n` +
    `bodies = ${bodies.length}\n` +
    `KE     = ${fmt(en.ke)}\n` +
    `PE     = ${fmt(en.pe)}\n` +
    `E      = ${fmt(en.e)}\n` +
    `drift  = ${drift.toFixed(4)} %`;

  drawEnergyGraph();
  drawCurvGraph();
  requestAnimationFrame(frame);
}

loadPreset("orbit");
updateOutputs();
requestAnimationFrame(frame);
