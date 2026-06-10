/* ============================================================
   Gravity Lab — physics core
   Internal units: AU, day, M☉  (so G = k², the Gauss constant²)
   Everything here is pure simulation logic, no DOM.
   ============================================================ */

"use strict";

// ---------------- units & constants ----------------

const U = {
  G_REAL: 2.959122082855911e-4,   // AU³ M☉⁻¹ day⁻²
  C_REAL: 173.14463267424,        // AU / day
  AU_KM: 1.495978707e8,
  MSUN_KG: 1.98892e30,
  KMS_PER_AUDAY: 1.495978707e8 / 86400,   // ≈ 1731.46 km/s
  YEAR_D: 365.25,
  ARCSEC: 206264.806,             // arcsec per radian
};

// ---------------- parameters ----------------

const P = {
  Gmult: 1,                 // multiplier on the real G
  cMult: 1,                 // multiplier on the real c
  n: 2,                     // force exponent (exploratory; GR terms assume 2)
  eps: 0,                   // Plummer softening (AU)
  dt: 0.5,                  // days
  sub: 8,                   // user-level steps per frame
  integrator: "verlet",     // euler | verlet | yoshida4 | rk4
  adaptive: true,           // subdivide steps on close encounters
  pn1: false,               // 1PN (EIH pairwise) relativistic correction
  pw: false,                // Paczyński–Wiita potential around the heaviest body
  merge: false,
};
const Gv = () => U.G_REAL * P.Gmult;
const Cv = () => U.C_REAL * P.cMult;

// ---------------- bodies ----------------

let bodies = [];
let time = 0;             // days
let bodyIdSeq = 1;

const PALETTE = ["#ffd166","#7c9bff","#ff6b8a","#5ce8c5","#c792ea","#ff9d5c","#7fd4ff","#aef07a","#f2f0a0","#8fa2ff"];
let colorIdx = 0;

function makeBody(m, x, y, vx, vy, color, name) {
  return {
    id: bodyIdSeq++,
    name: name || `body ${bodyIdSeq - 1}`,
    m, x, y, vx, vy,
    ax: 0, ay: 0,
    radAU: bodyRadiusAU(m),
    color: color || PALETTE[colorIdx++ % PALETTE.length],
    trail: [],
    // measurement state
    _pr1: null, _pr2: null, _pa1: 0, _pa2: 0,    // perihelion detector
    peri: [],                                    // [{t, angle}]
    _angAcc: 0, _angPrev: null, _wrapT: null,    // empirical period
    Tmeas: null,
  };
}

// stars/planets near solar density: R ≈ R☉·(m/M☉)^⅓; R☉ = 0.00465 AU
function bodyRadiusAU(m) { return 0.005 * Math.cbrt(Math.abs(m)); }

function heaviestIndex() {
  let k = -1, mm = -Infinity;
  for (let i = 0; i < bodies.length; i++) if (bodies[i].m > mm) { mm = bodies[i].m; k = i; }
  return k;
}
function heaviest() { return bodies[heaviestIndex()] || null; }

// Schwarzschild radius of a mass (AU)
function schwarzschildR(m) { const c = Cv(); return 2 * Gv() * m / (c * c); }

// ---------------- forces ----------------

// Newtonian (general n) + optional Paczyński–Wiita for pairs with the
// heaviest body + optional pairwise 1PN (EIH without third-body cross terms).
function computeAccels() {
  const N = bodies.length;
  for (const b of bodies) { b.ax = 0; b.ay = 0; }
  if (N < 2) return;
  const G = Gv(), c = Cv(), c2 = c * c;
  const e2 = P.eps * P.eps;
  const bh = P.pw ? heaviestIndex() : -1;
  const rsBH = bh >= 0 ? schwarzschildR(bodies[bh].m) : 0;

  for (let i = 0; i < N; i++) {
    const bi = bodies[i];
    for (let j = i + 1; j < N; j++) {
      const bj = bodies[j];
      const dx = bj.x - bi.x, dy = bj.y - bi.y;
      const r2 = dx * dx + dy * dy + e2;
      const r = Math.sqrt(r2);

      let f;   // acceleration prefactor: a_i = f · m_j · (dx, dy)
      if (bh === i || bh === j) {
        // PW: F = G m M / (r − r_s)², exact ISCO at 6GM/c²
        const rr = Math.max(r - rsBH, r * 1e-4);
        f = G / (rr * rr * r);
      } else if (P.n === 2) {
        f = G / (r2 * r);
      } else {
        f = G * Math.pow(r2, -(P.n + 1) / 2);
      }
      bi.ax += f * bj.m * dx;  bi.ay += f * bj.m * dy;
      bj.ax -= f * bi.m * dx;  bj.ay -= f * bi.m * dy;

      if (P.pn1) {
        // 1PN two-body (EIH) term. Convention: N̂ = (x_body − x_companion)/r,
        // i.e. pointing from the companion OUT to the accelerated body, so the
        // test-particle limit is Will's standard EOM
        //   a = (GM/c²r²)[ N̂(4GM/r − v²) + 4(N̂·v)v ]
        const inv = 1 / r;
        const gr2 = G / (c2 * r2);
        const vi2 = bi.vx * bi.vx + bi.vy * bi.vy;
        const vj2 = bj.vx * bj.vx + bj.vy * bj.vy;
        const vivj = bi.vx * bj.vx + bi.vy * bj.vy;
        // on i: N̂ = (x_i − x_j)/r = −(dx, dy)/r
        {
          const nx = -dx * inv, ny = -dy * inv;
          const nvi = nx * bi.vx + ny * bi.vy;
          const nvj = nx * bj.vx + ny * bj.vy;
          const A = 4 * G * bj.m * inv + 5 * G * bi.m * inv
                  - vi2 + 4 * vivj - 2 * vj2 + 1.5 * nvj * nvj;
          const B = 4 * nvi - 3 * nvj;
          bi.ax += gr2 * bj.m * (nx * A + (bi.vx - bj.vx) * B);
          bi.ay += gr2 * bj.m * (ny * A + (bi.vy - bj.vy) * B);
        }
        // on j: N̂ = (x_j − x_i)/r = +(dx, dy)/r, roles swapped
        {
          const nx = dx * inv, ny = dy * inv;
          const nvi = nx * bi.vx + ny * bi.vy;
          const nvj = nx * bj.vx + ny * bj.vy;
          const A = 4 * G * bi.m * inv + 5 * G * bj.m * inv
                  - vj2 + 4 * vivj - 2 * vi2 + 1.5 * nvi * nvi;
          const B = 4 * nvj - 3 * nvi;
          bj.ax += gr2 * bi.m * (nx * A + (bj.vx - bi.vx) * B);
          bj.ay += gr2 * bi.m * (ny * A + (bj.vy - bi.vy) * B);
        }
      }
    }
  }
}

// pair potential consistent with the softened force (Newtonian part)
function pairPotential(mi, mj, r, isBHpair, rsBH) {
  const G = Gv();
  if (isBHpair) return -G * mi * mj / Math.max(r - rsBH, r * 1e-4);
  if (Math.abs(P.n - 1) < 1e-9) return G * mi * mj * Math.log(r);
  return -G * mi * mj / ((P.n - 1) * Math.pow(r, P.n - 1));
}

// ---------------- conserved quantities ----------------

function energies() {
  let ke = 0, pe = 0;
  const e2 = P.eps * P.eps;
  const bh = P.pw ? heaviestIndex() : -1;
  const rsBH = bh >= 0 ? schwarzschildR(bodies[bh].m) : 0;
  for (let i = 0; i < bodies.length; i++) {
    const bi = bodies[i];
    ke += 0.5 * bi.m * (bi.vx * bi.vx + bi.vy * bi.vy);
    for (let j = i + 1; j < bodies.length; j++) {
      const bj = bodies[j];
      const dx = bj.x - bi.x, dy = bj.y - bi.y;
      pe += pairPotential(bi.m, bj.m, Math.sqrt(dx * dx + dy * dy + e2),
                          bh === i || bh === j, rsBH);
    }
  }
  return { ke, pe, e: ke + pe };
}

function momentum() {
  let px = 0, py = 0;
  for (const b of bodies) { px += b.m * b.vx; py += b.m * b.vy; }
  return [px, py];
}

function angularMomentum() {   // about the origin; conserved for isolated systems
  let L = 0;
  for (const b of bodies) L += b.m * (b.x * b.vy - b.y * b.vx);
  return L;
}

// gravitational potential per unit mass at a point (visuals + diagnostics);
// pass `skip` to exclude a body's own contribution
function potentialAt(x, y, skip) {
  const e2 = P.eps * P.eps;
  const G = Gv();
  const bh = P.pw ? heaviest() : null;
  const rsBH = bh ? schwarzschildR(bh.m) : 0;
  let phi = 0;
  for (const b of bodies) {
    if (b === skip) continue;
    const dx = b.x - x, dy = b.y - y;
    const r = Math.sqrt(dx * dx + dy * dy + e2);
    if (b === bh) phi += -G * b.m / Math.max(r - rsBH, r * 1e-4);
    else if (Math.abs(P.n - 1) < 1e-9) phi += G * b.m * Math.log(r);
    else phi += -G * b.m / ((P.n - 1) * Math.pow(r, P.n - 1));
  }
  return phi;
}

// weak-field clock rate dτ/dt, clamped at the horizon
function clockRate(phi, v2) {
  const c2 = Cv() * Cv();
  return Math.sqrt(Math.max(0, 1 + 2 * phi / c2 - v2 / c2));
}

// ---------------- integrators ----------------

function getState() {
  const s = new Float64Array(bodies.length * 4);
  bodies.forEach((b, i) => { s[i*4] = b.x; s[i*4+1] = b.y; s[i*4+2] = b.vx; s[i*4+3] = b.vy; });
  return s;
}
function setState(s) {
  bodies.forEach((b, i) => { b.x = s[i*4]; b.y = s[i*4+1]; b.vx = s[i*4+2]; b.vy = s[i*4+3]; });
}
function derivOf(s, out) {
  setState(s);
  computeAccels();
  bodies.forEach((b, i) => { out[i*4] = b.vx; out[i*4+1] = b.vy; out[i*4+2] = b.ax; out[i*4+3] = b.ay; });
}

const INTEGRATORS = {
  euler(h) {           // symplectic (semi-implicit) Euler — 1st order
    computeAccels();
    for (const b of bodies) {
      b.vx += b.ax * h; b.vy += b.ay * h;
      b.x  += b.vx * h; b.y  += b.vy * h;
    }
  },

  verlet(h) {          // velocity Verlet — 2nd order symplectic
    computeAccels();
    for (const b of bodies) {
      b.x += b.vx * h + 0.5 * b.ax * h * h;
      b.y += b.vy * h + 0.5 * b.ay * h * h;
      b._ox = b.ax; b._oy = b.ay;
    }
    if (P.pn1) {       // velocity-dependent force: one predictor pass
      for (const b of bodies) { b._sx = b.vx; b._sy = b.vy; b.vx += b._ox * h; b.vy += b._oy * h; }
      computeAccels();
      for (const b of bodies) { b.vx = b._sx; b.vy = b._sy; }
    } else {
      computeAccels();
    }
    for (const b of bodies) {
      b.vx += 0.5 * (b._ox + b.ax) * h;
      b.vy += 0.5 * (b._oy + b.ay) * h;
    }
  },

  yoshida4(h) {        // Yoshida 1990 — 4th order symplectic composition
    const cb = Math.cbrt(2);
    const w1 = 1 / (2 - cb), w0 = -cb * w1;
    const cs = [w1 / 2, (w0 + w1) / 2, (w0 + w1) / 2, w1 / 2];
    const ds = [w1, w0, w1];
    for (let k = 0; k < 4; k++) {
      for (const b of bodies) { b.x += cs[k] * b.vx * h; b.y += cs[k] * b.vy * h; }
      if (k < 3) {
        computeAccels();
        for (const b of bodies) { b.vx += ds[k] * b.ax * h; b.vy += ds[k] * b.ay * h; }
      }
    }
  },

  rk4(h) {             // classical Runge–Kutta — 4th order, handles v-dependent forces
    const n4 = bodies.length * 4;
    const s0 = getState();
    const k1 = new Float64Array(n4), k2 = new Float64Array(n4),
          k3 = new Float64Array(n4), k4 = new Float64Array(n4),
          tmp = new Float64Array(n4);
    derivOf(s0, k1);
    for (let i = 0; i < n4; i++) tmp[i] = s0[i] + 0.5 * h * k1[i];
    derivOf(tmp, k2);
    for (let i = 0; i < n4; i++) tmp[i] = s0[i] + 0.5 * h * k2[i];
    derivOf(tmp, k3);
    for (let i = 0; i < n4; i++) tmp[i] = s0[i] + h * k3[i];
    derivOf(tmp, k4);
    for (let i = 0; i < n4; i++) tmp[i] = s0[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    setState(tmp);
    computeAccels();
  },
};

// shortest dynamical timescale among all pairs (for adaptive stepping)
function minTimescale() {
  let tau = Infinity;
  const G = Gv();
  for (let i = 0; i < bodies.length; i++) {
    for (let j = i + 1; j < bodies.length; j++) {
      const a = bodies[i], b = bodies[j];
      const dx = b.x - a.x, dy = b.y - a.y;
      const r = Math.sqrt(dx * dx + dy * dy) + 1e-12;
      const dvx = b.vx - a.vx, dvy = b.vy - a.vy;
      const vrel = Math.sqrt(dvx * dvx + dvy * dvy);
      if (vrel > 0) tau = Math.min(tau, r / vrel);
      const mu = G * (a.m + b.m);
      if (mu > 0) tau = Math.min(tau, Math.sqrt(r * r * r / mu));
    }
  }
  return tau;
}

// advance one user-level step P.dt (subdivided if adaptive)
function advance(dt) {
  let steps = 1;
  if (P.adaptive && bodies.length > 1) {
    const tau = minTimescale();
    if (isFinite(tau)) steps = Math.min(512, Math.max(1, Math.ceil(dt / (0.08 * tau))));
  }
  P._steps = steps;
  const h = dt / steps;
  const integ = INTEGRATORS[P.integrator] || INTEGRATORS.verlet;
  for (let k = 0; k < steps; k++) {
    integ(h);
    time += h;
    trackOrbits();
  }
}

// ---------------- collisions & capture ----------------

let onBodyRemoved = null;   // hook for the UI (selection fixup)

function mergeCollisions() {
  for (let i = 0; i < bodies.length; i++) {
    for (let j = i + 1; j < bodies.length; j++) {
      const a = bodies[i], b = bodies[j];
      const dx = b.x - a.x, dy = b.y - a.y;
      const rc = a.radAU + b.radAU;
      if (dx * dx + dy * dy < rc * rc) {
        absorb(a, b);
        bodies.splice(j, 1); j--;
      }
    }
  }
}

// bodies (and photons) falling inside 2 r_s of a PW black hole are captured
function captureByBH() {
  const bh = heaviest();
  if (!bh) return;
  const rCap = 2 * schwarzschildR(bh.m);
  for (let j = bodies.length - 1; j >= 0; j--) {
    const b = bodies[j];
    if (b === bh) continue;
    const dx = b.x - bh.x, dy = b.y - bh.y;
    if (dx * dx + dy * dy < rCap * rCap) {
      absorb(bh, b);
      bodies.splice(j, 1);
    }
  }
}

function absorb(a, b) {     // a absorbs b, conserving mass & momentum
  const m = a.m + b.m;
  a.x = (a.x * a.m + b.x * b.m) / m;
  a.y = (a.y * a.m + b.y * b.m) / m;
  a.vx = (a.vx * a.m + b.vx * b.m) / m;
  a.vy = (a.vy * a.m + b.vy * b.m) / m;
  a.m = m;
  a.radAU = bodyRadiusAU(m);
  if (onBodyRemoved) onBodyRemoved(b, a);
}

// ---------------- orbit measurements ----------------
// Run at every substep: perihelion passages (precession) and empirical
// orbital periods (Kepler T²–a³ test), all relative to the heaviest body.

const keplerPoints = [];    // {id, name, color, a, T}

function trackOrbits() {
  const cIdx = heaviestIndex();
  if (cIdx < 0 || bodies.length < 2) return;
  const c = bodies[cIdx];
  for (const b of bodies) {
    if (b === c) continue;
    const rx = b.x - c.x, ry = b.y - c.y;
    const r = Math.sqrt(rx * rx + ry * ry);

    // perihelion: local minimum of r, refined to sub-step accuracy with a
    // parabolic fit so precession far smaller than Δθ per step is measurable
    const ang0 = Math.atan2(ry, rx);
    if (b._pr1 !== null && b._pr2 !== null && b._pr1 < b._pr2 && b._pr1 < r) {
      const denom = b._pr2 - 2 * b._pr1 + r;
      let off = denom > 0 ? 0.5 * (b._pr2 - r) / denom : 0;   // min offset in steps
      off = Math.max(-1, Math.min(1, off));
      let dAng = ang0 - b._pa2;
      while (dAng > Math.PI) dAng -= 2 * Math.PI;
      while (dAng < -Math.PI) dAng += 2 * Math.PI;
      b.peri.push({ t: time, angle: b._pa1 + off * dAng / 2 });
      if (b.peri.length > 60) b.peri.shift();
    }
    b._pr2 = b._pr1; b._pa2 = b._pa1;
    b._pr1 = r; b._pa1 = ang0;

    // empirical period: accumulate swept angle until 2π
    const ang = ang0;
    if (b._angPrev !== null) {
      let d = ang - b._angPrev;
      if (d > Math.PI) d -= 2 * Math.PI;
      if (d < -Math.PI) d += 2 * Math.PI;
      b._angAcc += d;
      if (Math.abs(b._angAcc) >= 2 * Math.PI) {
        if (b._wrapT !== null) {
          b.Tmeas = time - b._wrapT;
          const el = orbitalElements(b, c);
          if (el && el.a > 0) {
            const k = keplerPoints.findIndex(p => p.id === b.id);
            const pt = { id: b.id, name: b.name, color: b.color, a: el.a, T: b.Tmeas };
            if (k >= 0) keplerPoints[k] = pt; else keplerPoints.push(pt);
          }
        }
        b._wrapT = time;
        b._angAcc -= Math.sign(b._angAcc) * 2 * Math.PI;
      }
    } else {
      b._wrapT = time;
    }
    b._angPrev = ang;
  }
}

// osculating Keplerian elements of b about central mass c (Newtonian, n=2)
function orbitalElements(b, c) {
  if (!c || b === c) return null;
  const mu = Gv() * (b.m + c.m);
  const rx = b.x - c.x, ry = b.y - c.y;
  const vx = b.vx - c.vx, vy = b.vy - c.vy;
  const r = Math.sqrt(rx * rx + ry * ry);
  const v2 = vx * vx + vy * vy;
  const eOrb = v2 / 2 - mu / r;                  // specific orbital energy
  const h = rx * vy - ry * vx;                   // specific angular momentum
  const a = -mu / (2 * eOrb);
  const e = Math.sqrt(Math.max(0, 1 + 2 * eOrb * h * h / (mu * mu)));
  const T = a > 0 ? 2 * Math.PI * Math.sqrt(a * a * a / mu) : NaN;
  return { a, e, T, r, v: Math.sqrt(v2), mu, bound: eOrb < 0 };
}

// measured perihelion precession (rad/orbit, mean of recent passages)
function measuredPrecession(b) {
  if (!b || b.peri.length < 2) return null;
  const k = Math.max(1, b.peri.length - 6);
  let sum = 0, cnt = 0;
  for (let i = k; i < b.peri.length; i++) {
    let d = b.peri[i].angle - b.peri[i - 1].angle;
    while (d > Math.PI) d -= 2 * Math.PI;
    while (d < -Math.PI) d += 2 * Math.PI;
    sum += d; cnt++;
  }
  return cnt ? sum / cnt : null;
}

// 1PN prediction: Δϖ = 6πG(M+m) / (c² a (1−e²)) per orbit
function predictedPrecession(b, c) {
  const el = orbitalElements(b, c);
  if (!el || !el.bound) return null;
  const cc = Cv() * Cv();
  return 6 * Math.PI * el.mu / (cc * el.a * (1 - el.e * el.e));
}

// ---------------- photons (weak-field light rays) ----------------
// Massless rays at speed c; transverse acceleration is 2× Newtonian
// (the GR weak-field factor), |v| renormalized to c each step.

let photons = [];
let onPhotonDied = null;    // hook: UI keeps finished ray paths as fading ghosts

function emitPhotonFan(x0, y0, dirX, dirY, spread, count, width) {
  const c = Cv();
  const len = Math.hypot(dirX, dirY) || 1;
  const ux = dirX / len, uy = dirY / len;
  const px = -uy, py = ux;                       // perpendicular
  for (let i = 0; i < count; i++) {
    const off = (i / (count - 1) - 0.5) * width;
    photons.push({
      x: x0 + px * off, y: y0 + py * off,
      vx: ux * c, vy: uy * c,
      trail: [],
    });
  }
  void spread;
}

function stepPhotons(dt, killR2, sampleDist) {
  const c = Cv(), G = Gv();
  const e2 = P.eps * P.eps;
  const bh = P.pw ? heaviest() : null;
  const rCap = bh ? 1.5 * schwarzschildR(bh.m) : 0;
  const sd2 = sampleDist * sampleDist;
  for (let p = photons.length - 1; p >= 0; p--) {
    const ph = photons[p];
    if (!ph.trail.length) ph.trail.push(ph.x, ph.y);
    // micro-steps so the ray resolves the field near each mass; the trail is
    // sampled by distance travelled (a ray may cross the view in one frame)
    let remaining = dt;
    let dead = false;
    let guard = 4000;
    while (remaining > 0 && !dead && guard-- > 0) {
      let rMin = Infinity, ax = 0, ay = 0;
      for (const b of bodies) {
        const dx = b.x - ph.x, dy = b.y - ph.y;
        const r2 = dx * dx + dy * dy + e2;
        const r = Math.sqrt(r2);
        rMin = Math.min(rMin, r);
        const f = 2 * G * b.m / (r2 * r);        // GR factor 2
        ax += f * dx; ay += f * dy;
      }
      const h = Math.min(remaining, Math.max(1e-9, 0.05 * rMin / c));
      ph.vx += ax * h; ph.vy += ay * h;
      const s = c / Math.hypot(ph.vx, ph.vy);    // keep |v| = c
      ph.vx *= s; ph.vy *= s;
      ph.x += ph.vx * h; ph.y += ph.vy * h;
      remaining -= h;
      const n = ph.trail.length;
      const lx = ph.x - ph.trail[n - 2], ly = ph.y - ph.trail[n - 1];
      if (lx * lx + ly * ly > sd2) {
        ph.trail.push(ph.x, ph.y);
        if (ph.trail.length > 1200) ph.trail.splice(0, ph.trail.length - 1200);
      }
      if (bh) {
        const dx = ph.x - bh.x, dy = ph.y - bh.y;
        if (dx * dx + dy * dy < rCap * rCap) dead = true;   // photon capture
      }
      if (ph.x * ph.x + ph.y * ph.y > killR2) break;
    }
    if (dead || ph.x * ph.x + ph.y * ph.y > killR2) {
      ph.trail.push(ph.x, ph.y);
      if (onPhotonDied) onPhotonDied(ph.trail);
      photons.splice(p, 1);
    }
  }
}

// ---------------- Lagrange points (circular restricted 3-body) ----------------
// Returns L1..L5 in *inertial* coordinates for the two heaviest bodies,
// plus the frame info {bx, by, theta, omega, d}.

function lagrangeData() {
  if (bodies.length < 2) return null;
  const sorted = [...bodies].sort((a, b) => b.m - a.m);
  const b1 = sorted[0], b2 = sorted[1];
  const dx = b2.x - b1.x, dy = b2.y - b1.y;
  const d = Math.hypot(dx, dy);
  if (d <= 0) return null;
  const G = Gv();
  const M = b1.m + b2.m;
  const omega = Math.sqrt(G * M / (d * d * d));
  const theta = Math.atan2(dy, dx);
  const bx = (b1.x * b1.m + b2.x * b2.m) / M;
  const by = (b1.y * b1.m + b2.y * b2.m) / M;
  // frame coords: barycenter origin, b2 on +x
  const x1 = -d * b2.m / M, x2 = d * b1.m / M;
  const w2 = omega * omega;
  const fx = x => {     // net axial force in the rotating frame (y = 0)
    const r1 = x - x1, r2_ = x - x2;
    return -G * b1.m * Math.sign(r1) / (r1 * r1)
           - G * b2.m * Math.sign(r2_) / (r2_ * r2_)
           + w2 * x;
  };
  const bisect = (lo, hi) => {
    let flo = fx(lo);
    for (let it = 0; it < 80; it++) {
      const mid = (lo + hi) / 2, fm = fx(mid);
      if (flo * fm <= 0) hi = mid; else { lo = mid; flo = fm; }
    }
    return (lo + hi) / 2;
  };
  const eps = d * 1e-6;
  const L = [
    bisect(x1 + eps, x2 - eps),            // L1
    bisect(x2 + eps, x2 + 3 * d),          // L2
    bisect(x1 - 3 * d, x1 - eps),          // L3
  ].map(x => [x, 0]);
  L.push([x1 + d / 2, +Math.sqrt(3) / 2 * d]);   // L4
  L.push([x1 + d / 2, -Math.sqrt(3) / 2 * d]);   // L5
  const ct = Math.cos(theta), st = Math.sin(theta);
  const inertial = L.map(([x, y]) => [bx + x * ct - y * st, by + x * st + y * ct]);
  return { points: inertial, frame: { bx, by, theta, omega, d, b1, b2 } };
}

// effective potential in the rotating frame (for grid coloring / profile)
function effectivePotentialAt(x, y, frame) {
  const rho2 = (x - frame.bx) ** 2 + (y - frame.by) ** 2;
  return potentialAt(x, y) - 0.5 * frame.omega * frame.omega * rho2;
}
