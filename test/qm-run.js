#!/usr/bin/env node
/* ============================================================
   Quantum Lab validation suite — run with:  node test/qm-run.js
   Loads the actual TDSE engine + scenarios headlessly and checks
   them against analytic quantum mechanics.
   ============================================================ */

"use strict";

const fs = require("fs");
const path = require("path");
const vm = require("vm");

const sandbox = { console, Math, performance: { now: () => 0 } };
vm.createContext(sandbox);
for (const f of ["qm-engine.js", "qm-scenarios.js"]) {
  vm.runInContext(fs.readFileSync(path.join(__dirname, "..", f), "utf8"), sandbox, { filename: f });
}
function S(code) { return vm.runInContext(`{${code}}`, sandbox); }

let failures = 0;
function check(name, cond, detail) {
  console.log(`${cond ? "PASS" : "FAIL"}  ${name}${detail ? `  [${detail}]` : ""}`);
  if (!cond) failures++;
}
function approx(a, b, tol) { return Math.abs(a - b) <= tol * Math.abs(b); }

// ---------- 1. FFT round-trip ----------
{
  const err = S(`
    const n = 256;
    const re = new Float64Array(n), im = new Float64Array(n);
    for (let i = 0; i < n; i++) { re[i] = Math.sin(i * 0.7) + 0.3 * i / n; im[i] = Math.cos(i * 1.3); }
    const r0 = re.slice(), i0 = im.slice();
    fft(re, im, false); fft(re, im, true);
    let e = 0;
    for (let i = 0; i < n; i++) e = Math.max(e, Math.abs(re[i] - r0[i]), Math.abs(im[i] - i0[i]));
    e;
  `);
  check("FFT round-trip exact to machine precision", err < 1e-12, `max err = ${err.toExponential(1)}`);
}

// ---------- 2. unitarity & energy conservation ----------
{
  const r = S(`
    const s = new QMState({ dim: 1, n: 1024, L: 120 });
    s.setPotential(x => 0.125 * x * x);
    s.gaussian({ x0: 8, sx: 2, kx: 0.5 });
    s.normalize();
    const E0 = s.energy().E;
    s.step(0.01, 5000);
    ({ dN: Math.abs(s.norm() - 1), dE: Math.abs((s.energy().E - E0) / E0) });
  `);
  check("norm conserved to 1e-12 over 5000 steps", r.dN < 1e-12, `|ΔN| = ${r.dN.toExponential(1)}`);
  check("⟨H⟩ bounded to 1e-6 (split-step O(dt²), non-secular)", r.dE < 1e-6, `|ΔE/E| = ${r.dE.toExponential(1)}`);
}

// ---------- 3. free-packet dispersion ----------
{
  const r = S(`
    const s = new QMState({ dim: 1, n: 2048, L: 800 });
    s.setPotential(() => 0);
    const s0 = 5;
    s.gaussian({ x0: -150, sx: s0, kx: 1 });
    s.normalize();
    s.step(0.05, 4000);                         // t = 200
    const t = s.t;
    ({ meas: Math.sqrt(s.moments().varX),
       pred: s0 * Math.sqrt(1 + (t / (2 * s0 * s0)) ** 2),
       x: s.moments().x, xPred: -150 + 1 * t });
  `);
  check("free packet spreads as σ₀√(1+(t/2σ₀²)²)", approx(r.meas, r.pred, 0.005),
        `σ = ${r.meas.toFixed(3)} vs ${r.pred.toFixed(3)}`);
  check("⟨x⟩ moves at group velocity k₀", approx(r.x, r.xPred, 0.005),
        `⟨x⟩ = ${r.x.toFixed(2)} vs ${r.xPred.toFixed(2)}`);
}

// ---------- 4. harmonic oscillator ----------
{
  const r = S(`
    const s = new QMState({ dim: 1, n: 1024, L: 120 });
    const w = 0.5;
    s.setPotential(x => 0.5 * w * w * x * x);
    const eig = solveEigen(s, 3, { dtau: 2e-3 });
    // coherent state: full revival after one classical period
    s.gaussian({ x0: 12, sx: Math.sqrt(0.5 / w) });
    s.normalize();
    const r0 = s.re.slice(), i0 = s.im.slice();
    const T = 2 * Math.PI / w;
    const steps = 4000;
    s.step(T / steps, steps);
    let ovR = 0, ovI = 0;
    for (let i = 0; i < s.n; i++) {
      ovR += r0[i] * s.re[i] + i0[i] * s.im[i];
      ovI += r0[i] * s.im[i] - i0[i] * s.re[i];
    }
    ({ E: eig.energies, w,
       overlap: Math.hypot(ovR, ovI) * s.dx,
       sig: Math.sqrt(s.moments().varX), sig0: Math.sqrt(0.5 / w) });
  `);
  check("HO eigenenergies = (n+½)ω via imaginary time",
        approx(r.E[0], 0.25, 0.002) && approx(r.E[1], 0.75, 0.002) && approx(r.E[2], 1.25, 0.002),
        r.E.map(e => e.toFixed(4)).join(", "));
  check("coherent state revives after T = 2π/ω", r.overlap > 0.9999,
        `|⟨ψ(0)|ψ(T)⟩| = ${r.overlap.toFixed(6)}`);
  check("coherent state does not spread", approx(r.sig, r.sig0, 0.01),
        `σ = ${r.sig.toFixed(4)} vs ${r.sig0.toFixed(4)}`);
}

// ---------- 5. quantum bouncer = Airy zeros (real qBOUNCE energies) ----------
{
  const r = S(`
    const s = new QMState({ dim: 1, n: 1024, L: 60 });
    s.setPotential(x => Math.abs(x));
    s.oddProject = true;
    const eig = solveEigen(s, 4, { dtau: 2e-3 });
    eig.energies;
  `);
  const scale = Math.cbrt(0.5);
  const exact = [2.33810741, 4.08794944, 5.52055983, 6.78670809].map(a => a * scale);
  const ok = r.every((E, n) => approx(E, exact[n], 0.005));
  check("bouncer eigenvalues match Airy zeros (E₁ ↔ 1.41 peV for neutrons)", ok,
        r.map((E, n) => `${E.toFixed(4)}/${exact[n].toFixed(4)}`).join(" "));
}

// ---------- 6. tunneling vs analytic transmission ----------
{
  const r = S(`
    const sc = QM_SCENARIOS.tunnel;
    const p = { V0: 1.0, a: 2, k0: 1.2 };
    const s = new QMState({ dim: 1, n: 2048, L: 600 });
    sc.build(s, p);
    s.step(0.04, 3750);                          // t = 150: cleared, before wrap-around
    ({ meas: s.probIn(p.a / 2 + 5, 1e9), pred: packetT(p.k0, 10, p.V0, sc.effWidth(s, p)) });
  `);
  check("transmission matches ∫|φ(k)|²T(k)dk within 2%", approx(r.meas, r.pred, 0.02),
        `T = ${r.meas.toFixed(4)} vs ${r.pred.toFixed(4)}`);
}

// ---------- 7. gravimeter phase = mgΔh·t/ħ ----------
{
  const r = S(`
    const sc = QM_SCENARIOS.gravimeter;
    const p = { g: 0.05, dh: 20 };
    const s = new QMState({ dim: 1, n: 1024, L: 120 });
    sc.build(s, p);
    for (let i = 0; i < 250; i++) { s.step(0.01, 8); sc.tick(s, p); }
    ({ slope: sc._acc / s.t, pred: p.g * p.dh });
  `);
  check("interferometer phase rate = mgΔh/ħ within 1%", approx(r.slope, r.pred, 0.01),
        `${r.slope.toFixed(5)} vs ${r.pred.toFixed(5)} rad/t`);
}

// ---------- 8. double slit fringes = λD/d ----------
{
  const r = S(`
    const sc = QM_SCENARIOS.doubleslit;
    const p = { k0: 2, d: 10, sw: 3 };
    const s = new QMState({ dim: 2, nx: 256, ny: 256, Lx: 120, Ly: 120 });
    sc.build(s, p);
    for (let i = 0; i < 1000; i++) { s.step(0.02, 4); sc.tick(s, p); }   // t = 80
    ({ meas: sc.fringeSpacing(s), pred: (2 * Math.PI / p.k0) * sc.screenX / p.d });
  `);
  check("first-minima fringe spacing matches λD/d within 5%", approx(r.meas, r.pred, 0.05),
        `${r.meas.toFixed(2)} vs ${r.pred.toFixed(2)}`);
}

// ---------- 9. Ehrenfest exact for quadratic potentials ----------
{
  const r = S(`
    const s = new QMState({ dim: 2, nx: 256, ny: 256, Lx: 80, Ly: 80 });
    const w = 0.3;
    s.setPotential((x, y) => 0.5 * w * w * (x * x + y * y));
    s.gaussian({ x0: 10, y0: 0, sx: 3, sy: 3, kx: 0, ky: 10 * w });   // circular classical orbit
    s.normalize();
    const T = 2 * Math.PI / w;
    s.step(T / 2000, 2000);
    const m = s.moments();
    // classical: returns to (10, 0) after one period
    ({ ex: m.x, ey: m.y });
  `);
  check("Ehrenfest: ⟨r⟩ follows the classical orbit exactly (harmonic)",
        Math.hypot(r.ex - 10, r.ey) < 0.05, `returned to (${r.ex.toFixed(3)}, ${r.ey.toFixed(3)})`);
}

// ---------- 10. entanglement from scattering ----------
{
  const r = S(`
    const sc = QM_SCENARIOS.entangle;
    const p = { V0: 2.5, k: 1.58 };              // E_rel = k² = V₀ → 50/50 split
    const s = new QMState({ dim: 2, nx: 128, ny: 128, Lx: 100, Ly: 100 });
    sc.build(s, p);
    const K0 = 1 / s.purity();
    s.step(0.02, 1600);                          // t = 32: well past the collision
    ({ K0, K1: 1 / s.purity() });
  `);
  check("initial product state has Schmidt number K = 1", Math.abs(r.K0 - 1) < 0.01,
        `K₀ = ${r.K0.toFixed(4)}`);
  check("50/50 collision gives maximal Schmidt number K ≈ 2", Math.abs(r.K1 - 2) < 0.15,
        `K after collision = ${r.K1.toFixed(3)} (two Schmidt modes: bounced / passed)`);
}

// ---------- 11. measurement-free no-interaction control ----------
{
  const r = S(`
    const sc = QM_SCENARIOS.entangle;
    const p = { V0: 0, k: 1.58 };
    const s = new QMState({ dim: 2, nx: 128, ny: 128, Lx: 100, Ly: 100 });
    sc.build(s, p);
    s.step(0.02, 1600);
    1 / s.purity();
  `);
  check("without interaction the particles stay unentangled (K ≈ 1)",
        Math.abs(r - 1) < 0.02, `K = ${r.toFixed(4)}`);
}

console.log(failures ? `\n${failures} FAILURE(S)` : "\nALL CHECKS PASSED");
process.exit(failures ? 1 : 0);
