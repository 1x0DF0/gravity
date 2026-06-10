#!/usr/bin/env node
/* ============================================================
   Gravity Lab validation suite — run with:  node test/run.js
   Loads the actual simulation code (physics.js + presets.js)
   headlessly and checks it against analytic physics.
   ============================================================ */

"use strict";

const fs = require("fs");
const path = require("path");
const vm = require("vm");

// minimal globals the sim files expect (presets touch V at build time)
const sandbox = {
  console,
  performance: { now: () => 0 },
  V: { scale: 1, rotFrame: false },
};
vm.createContext(sandbox);
for (const f of ["physics.js", "presets.js"]) {
  vm.runInContext(fs.readFileSync(path.join(__dirname, "..", f), "utf8"), sandbox, { filename: f });
}

// run test code inside the same context so it sees the sim's globals
// (wrapped in a block so const/let don't leak between snippets)
function S(code) { return vm.runInContext(`{${code}}`, sandbox); }

let failures = 0;
function check(name, cond, detail) {
  const ok = !!cond;
  console.log(`${ok ? "PASS" : "FAIL"}  ${name}${detail ? `  [${detail}]` : ""}`);
  if (!ok) failures++;
}
function approx(a, b, tol) { return Math.abs(a - b) <= tol * Math.abs(b); }

// preset loader mirroring main.js's reset semantics
S(`
globalThis.loadP = function (name) {
  bodies = []; photons = []; colorIdx = 0; time = 0; bodyIdSeq = 1;
  keplerPoints.length = 0;
  Object.assign(P, { Gmult: 1, cMult: 1, n: 2, eps: 0, dt: 0.5, sub: 12,
    integrator: "verlet", adaptive: true, pn1: false, pw: false, merge: false });
  PRESETS[name].build();
  computeAccels();
};
globalThis.run = function (days) {
  const steps = Math.ceil(days / P.dt);
  for (let i = 0; i < steps; i++) advance(P.dt);
};
`);

// ---------- 1. units ----------
{
  const G = S("U.G_REAL");
  check("G in AU³ M☉⁻¹ day⁻² equals the Gauss constant²",
        approx(G, 0.01720209895 ** 2, 1e-9), `G = ${G}`);
  const vEarth = S("Math.sqrt(Gv() * 1 / 1) * U.KMS_PER_AUDAY");
  check("circular speed at 1 AU around 1 M☉ ≈ 29.78 km/s",
        approx(vEarth, 29.78, 0.005), `${vEarth.toFixed(2)} km/s`);
  const cKms = S("U.C_REAL * U.KMS_PER_AUDAY");
  check("c in AU/day converts to 299792 km/s",
        approx(cKms, 299792.458, 1e-4), `${cKms.toFixed(0)} km/s`);
}

// ---------- 2. solar system sanity ----------
{
  S("loadP('solar')");
  const el = S(`
    const earth = bodies.find(b => b.name === "Earth");
    orbitalElements(earth, bodies.find(b => b.name === "Sun"));
  `);
  check("Earth osculating a ≈ 1 AU", approx(el.a, 1.0, 0.01), `a = ${el.a.toFixed(4)} AU`);
  check("Earth osculating period ≈ 365.25 d", approx(el.T, 365.25, 0.01), `T = ${el.T.toFixed(2)} d`);
  check("Earth eccentricity ≈ 0.0167", Math.abs(el.e - 0.0167) < 0.002, `e = ${el.e.toFixed(4)}`);
}

// ---------- 3. conservation (symplectic integrators) ----------
{
  // symplectic ⇒ bounded (non-secular) energy error; tighter for higher order
  for (const [integ, tol] of [["verlet", 1e-5], ["yoshida4", 1e-8]]) {
    S(`loadP('solar'); P.integrator = '${integ}'`);
    const r = S(`
      const e0 = energies().e, L0 = angularMomentum();
      run(2000);
      ({ dE: Math.abs((energies().e - e0) / e0),
         dL: Math.abs((angularMomentum() - L0) / L0) });
    `);
    check(`${integ}: energy drift < ${tol} over 2000 d`, r.dE < tol, `ΔE/E = ${r.dE.toExponential(2)}`);
    check(`${integ}: ang. momentum drift < 1e-9`, r.dL < 1e-9, `ΔL/L = ${r.dL.toExponential(2)}`);
  }
  // order comparison: yoshida4 should beat verlet at identical dt
  const drift = i => { S(`loadP('binary'); P.integrator='${i}'; P.adaptive=false; P.dt=0.5`);
    return S("const e0=energies().e; run(3000); Math.abs((energies().e-e0)/e0)"); };
  const dv = drift("verlet"), dy = drift("yoshida4");
  check("yoshida4 error ≪ verlet at same Δt", dy < dv / 30,
        `verlet ${dv.toExponential(2)} vs yoshida4 ${dy.toExponential(2)}`);
}

// ---------- 4. figure-8 choreography (real units) ----------
{
  S("loadP('figure8'); P.integrator='yoshida4'; P.adaptive=false; P.dt=0.02");
  const err = S(`
    const T = 6.32591398 / Math.sqrt(Gv());     // period rescaled to days
    const start = bodies.map(b => [b.x, b.y]);
    run(T);
    Math.max(...bodies.map((b, i) => Math.hypot(b.x - start[i][0], b.y - start[i][1])));
  `);
  check("figure-8 returns to start after one period", err < 0.02, `max err = ${err.toFixed(5)} AU`);
}

// ---------- 5. 1PN perihelion precession vs analytic ----------
{
  S("loadP('mercury'); P.dt = 0.1");
  const r = S(`
    run(88 * 12);                                // ~12 Mercury orbits
    const m = bodies.find(b => b.name === "Mercury");
    ({ meas: measuredPrecession(m), pred: predictedPrecession(m, bodies[0]), n: m.peri.length });
  `);
  check("Mercury 1PN precession matches 6πGM/(c²a(1−e²)) within 3%",
        r.meas !== null && approx(r.meas, r.pred, 0.03),
        `measured ${(r.meas * 206264.8).toFixed(1)}″ vs predicted ${(r.pred * 206264.8).toFixed(1)}″ per orbit`);
  // and with GR off it should vanish
  S("loadP('mercury'); P.pn1 = false; P.dt = 0.1");
  const r0 = S(`
    run(88 * 12);
    measuredPrecession(bodies.find(b => b.name === "Mercury"));
  `);
  check("precession ≈ 0 with 1PN off", Math.abs(r0) < Math.abs(r.pred) / 20,
        `${(r0 * 206264.8).toFixed(3)}″/orbit`);
}

// ---------- 6. Paczyński–Wiita ISCO ----------
{
  const orbitR = r => S(`
    loadP('isco');
    bodies = bodies.slice(0, 1);                 // keep only the BH
    const bh = bodies[0];
    const rg = Gv() * bh.m / (Cv() * Cv());
    const r0 = ${r} * rg;
    const rs = schwarzschildR(bh.m);
    const v = Math.sqrt(Gv() * bh.m * r0) / (r0 - rs);
    bodies.push(makeBody(1e-12, r0, 0, 0, v, "#fff", "probe"));
    computeAccels();
    P.integrator = "rk4"; P.dt = 1e-4;
    const T = 2 * Math.PI * r0 / v;
    run(5 * T);
    Math.hypot(bodies[1].x, bodies[1].y) / r0;
  `);
  const stable = orbitR(8), plunge = orbitR(5.2);
  check("PW circular orbit at 8 GM/c² stays circular (outside ISCO)",
        Math.abs(stable - 1) < 0.05, `r/r₀ = ${stable.toFixed(3)} after 5 orbits`);
  check("PW circular orbit at 5.2 GM/c² is unstable (inside ISCO at 6)",
        plunge < 0.5 || plunge > 2 || !isFinite(plunge), `r/r₀ = ${plunge.toFixed(3)}`);
}

// ---------- 7. light deflection ----------
{
  const r = S(`
    bodies = []; photons = []; time = 0;
    Object.assign(P, { Gmult: 1, cMult: 0.05, n: 2, eps: 0, pn1: false, pw: false });
    bodies.push(makeBody(1, 0, 0, 0, 0, "#fff", "Sun"));
    const b = 0.05;                              // impact parameter (AU)
    emitPhotonFan(-6, b, 1, 0, 0, 2, 0);         // two coincident rays @ y = b
    for (let i = 0; i < 400 && photons.length && photons[0].x < 6; i++) {
      stepPhotons(0.05, 1e9, 0.05);
    }
    const ph = photons[0];
    ({ defl: Math.atan2(-ph.vy, ph.vx), pred: 4 * Gv() / (Cv() * Cv() * b) });
  `);
  check("photon deflection matches 4GM/(c²b) within 5%",
        approx(r.defl, r.pred, 0.05),
        `measured ${(r.defl * 206264.8).toFixed(0)}″ vs predicted ${(r.pred * 206264.8).toFixed(0)}″`);
}

// ---------- 8. Lagrange points ----------
{
  S("loadP('trojans')");
  const r = S(`
    const d = lagrangeData();
    const jup = bodies.find(b => b.name === "Jupiter");
    const sun = bodies.find(b => b.name === "Sun");
    const [l1, , , l4] = [d.points[0], d.points[1], d.points[2], d.points[3]];
    ({ l1j: Math.hypot(l1[0] - jup.x, l1[1] - jup.y),
       l4s: Math.hypot(l4[0] - sun.x, l4[1] - sun.y),
       l4j: Math.hypot(l4[0] - jup.x, l4[1] - jup.y),
       aJ: Math.hypot(jup.x - sun.x, jup.y - sun.y) });
  `);
  const hill = r.aJ * Math.cbrt(9.546e-4 / 3);
  check("L1 sits ≈ one Hill radius from Jupiter", approx(r.l1j, hill, 0.10),
        `${r.l1j.toFixed(3)} AU vs Hill ${hill.toFixed(3)} AU`);
  check("L4 is equilateral (|L4−Sun| ≈ |L4−Jup| ≈ a)",
        approx(r.l4s, r.aJ, 0.01) && approx(r.l4j, r.aJ, 0.01),
        `${r.l4s.toFixed(3)} / ${r.l4j.toFixed(3)} vs ${r.aJ.toFixed(3)}`);
  // a probe started exactly at L4 with co-rotation velocity should stay near L4
  const drift = S(`
    const d = lagrangeData();
    const [x, y] = d.points[3];
    const f = d.frame;
    bodies.push(makeBody(1e-12, x, y,
      -f.omega * (y - f.by), f.omega * (x - f.bx), "#fff", "L4probe"));
    computeAccels();
    run(4332 * 2);                               // two Jupiter periods
    const d2 = lagrangeData();
    const p = bodies.find(b => b.name === "L4probe");
    Math.hypot(p.x - d2.points[3][0], p.y - d2.points[3][1]);
  `);
  check("probe parked at L4 librates (stays within 0.5 AU over 2 Jupiter orbits)",
        drift < 0.5, `drift = ${drift.toFixed(3)} AU`);
}

// ---------- 9. time dilation ----------
{
  S("loadP('solar')");
  const r = S(`
    const e = bodies.find(b => b.name === "Earth");
    const phi = potentialAt(e.x, e.y, e);
    const v2 = e.vx * e.vx + e.vy * e.vy;
    1 - clockRate(phi, v2);
  `);
  // GM/(r c²) + v²/2c² ≈ 9.87e-9 + 4.93e-9 ≈ 1.48e-8 at Earth
  check("Earth clock runs slow by ≈ 1.48×10⁻⁸", approx(r, 1.48e-8, 0.05),
        `1−dτ/dt = ${r.toExponential(3)}`);
}

// ---------- 10. merging conserves momentum ----------
{
  const r = S(`
    bodies = []; photons = []; time = 0;
    Object.assign(P, { Gmult: 1, cMult: 1, n: 2, eps: 0, pn1: false, pw: false });
    bodies.push(makeBody(2e-6, 0, 0, 0.001, 0.002, "#fff", "a"));
    bodies.push(makeBody(3e-6, 1e-6, 0, -0.003, 0.001, "#fff", "b"));
    const p0 = momentum();
    mergeCollisions();
    const p1 = momentum();
    ({ n: bodies.length, err: Math.hypot(p1[0] - p0[0], p1[1] - p0[1]), m: bodies[0].m });
  `);
  check("contact merge conserves momentum & mass", r.n === 1 && r.err < 1e-20 && approx(r.m, 5e-6, 1e-12),
        `|Δp| = ${r.err.toExponential(1)}`);
}

// ---------- 11. empirical Kepler periods ----------
{
  S("loadP('solar'); P.integrator='yoshida4'");
  const r = S(`
    run(700);                                    // > 1 orbit for the inner planets
    const pts = keplerPoints.filter(p => ["Mercury","Venus","Earth"].includes(p.name));
    pts.map(p => ({ name: p.name, ratio: p.T / (2 * Math.PI * Math.sqrt(p.a ** 3 / Gv())) }));
  `);
  check("measured periods of Mercury/Venus/Earth match 2π√(a³/GM) within 1%",
        r.length === 3 && r.every(p => Math.abs(p.ratio - 1) < 0.01),
        r.map(p => `${p.name} ${p.ratio.toFixed(4)}`).join(", "));
}

console.log(failures ? `\n${failures} FAILURE(S)` : "\nALL CHECKS PASSED");
process.exit(failures ? 1 : 0);
