/* ============================================================
   Gravity Lab — scenarios
   Real units throughout: AU, day, M☉.
   Solar-system bodies use real masses, semi-major axes and
   eccentricities, started at perihelion.
   ============================================================ */

"use strict";

// ---------- seeded RNG (reproducibility) ----------
let rngSeed = 42;
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// place a body at perihelion of an orbit (a, e) about `central`, rotated by ang
function atPerihelion(central, m, a, e, ang, color, name) {
  const mu = Gv() * (central.m + m);
  const rp = a * (1 - e);
  const vp = Math.sqrt(mu * (1 + e) / rp);
  const ca = Math.cos(ang), sa = Math.sin(ang);
  return makeBody(m,
    central.x + rp * ca, central.y + rp * sa,
    central.vx - vp * sa, central.vy + vp * ca,
    color, name);
}

// real planetary data: [name, a (AU), e, mass (M☉), color]
const PLANETS = [
  ["Mercury", 0.3871,  0.2056, 1.660e-7, "#b9a48a"],
  ["Venus",   0.7233,  0.0068, 2.448e-6, "#e8c97a"],
  ["Earth",   1.0000,  0.0167, 3.003e-6, "#6fa8ff"],
  ["Mars",    1.5237,  0.0934, 3.227e-7, "#ff7a5c"],
  ["Jupiter", 5.2026,  0.0489, 9.546e-4, "#e0b080"],
  ["Saturn",  9.5549,  0.0565, 2.858e-4, "#e8d9a0"],
  ["Uranus", 19.218,   0.0460, 4.366e-5, "#9fdcdc"],
  ["Neptune",30.110,   0.0095, 5.151e-5, "#7c9bff"],
];

const PRESETS = {

  solar: {
    label: "Solar system (real data)",
    info: "Real masses, semi-major axes and eccentricities; bodies start at perihelion.",
    build() {
      P.dt = 0.5; P.sub = 12; P.eps = 0;
      P.integrator = "yoshida4"; P.adaptive = true;
      V.scale = 160;
      const sun = makeBody(1, 0, 0, 0, 0, "#ffd166", "Sun");
      bodies.push(sun);
      const rnd = mulberry32(rngSeed);
      for (const [name, a, e, m, col] of PLANETS) {
        bodies.push(atPerihelion(sun, m, a, e, rnd() * 2 * Math.PI, col, name));
      }
      zeroTotalMomentum();
    },
  },

  mercury: {
    label: "Mercury precession (GR test)",
    info: "Sun + Mercury with the 1PN correction. c is scaled ×0.02 so the relativistic " +
          "perihelion advance (×2500 the real 0.104″/orbit) is visible and measurable. " +
          "Set the c multiplier to 1 for the real value.",
    build() {
      P.dt = 0.25; P.sub = 16; P.eps = 0;
      P.integrator = "rk4"; P.adaptive = true;
      P.pn1 = true; P.cMult = 0.02;
      V.scale = 700;
      const sun = makeBody(1, 0, 0, 0, 0, "#ffd166", "Sun");
      bodies.push(sun);
      bodies.push(atPerihelion(sun, 1.660e-7, 0.3871, 0.2056, 0, "#b9a48a", "Mercury"));
      zeroTotalMomentum();
    },
    select: "Mercury",
  },

  s2: {
    label: "S2 around Sgr A* (real orbit)",
    info: "The star S2 on its real 16-year, e = 0.88 orbit around the 4.3-million-solar-mass " +
          "black hole at the Galactic Centre, with the 1PN correction and real c. " +
          "GRAVITY measured its Schwarzschild precession: ≈ 12′ per orbit.",
    build() {
      P.dt = 2; P.sub = 60; P.eps = 0;
      P.integrator = "rk4"; P.adaptive = true;
      P.pn1 = true; P.cMult = 1;
      V.scale = 0.4;
      const bh = makeBody(4.297e6, 0, 0, 0, 0, "#12121f", "Sgr A*");
      bodies.push(bh);
      bodies.push(atPerihelion(bh, 13.6, 1031, 0.884, 1.0, "#7fd4ff", "S2"));
      bodies.push(atPerihelion(bh, 10,   1500, 0.70, 3.5, "#c792ea", "S-star b"));
      zeroTotalMomentum();
    },
    select: "S2",
  },

  isco: {
    label: "Black hole ISCO (Paczyński–Wiita)",
    info: "Test stars around Sgr A* in the Paczyński–Wiita potential Φ = −GM/(r−rₛ), " +
          "which reproduces the Schwarzschild ISCO at 6GM/c² (≈ 0.25 AU here). " +
          "The inner orbits are unstable and plunge; bodies inside 2rₛ are captured.",
    build() {
      P.dt = 2e-3; P.sub = 24; P.eps = 0;
      P.integrator = "rk4"; P.adaptive = true;
      P.pw = true; P.cMult = 1;
      V.scale = 250;
      const bh = makeBody(4.297e6, 0, 0, 0, 0, "#12121f", "Sgr A*");
      bodies.push(bh);
      const G = Gv(), rs = schwarzschildR(bh.m);
      const ring = [[1.0, "#7fd4ff"], [0.6, "#5ce8c5"], [0.35, "#ffd166"], [0.27, "#ff6b8a"]];
      ring.forEach(([r, col], i) => {
        const v = Math.sqrt(G * bh.m * r) / (r - rs);   // PW circular speed
        const a = i * 1.7;
        bodies.push(makeBody(1e-9, r * Math.cos(a), r * Math.sin(a),
          -v * Math.sin(a), v * Math.cos(a), col, `star r=${r} AU`));
      });
      // one eccentric plunger
      bodies.push(atPerihelion(bh, 1e-9, 0.8, 0.85, 4.0, "#ff9d5c", "plunger"));
    },
  },

  trojans: {
    label: "Sun–Jupiter Trojans (L4/L5)",
    info: "Test asteroids near Jupiter's L4 and L5 points. Switch on the rotating frame " +
          "to watch them librate around the Lagrange points (tadpole orbits).",
    build() {
      P.dt = 2; P.sub = 16; P.eps = 0;
      P.integrator = "yoshida4"; P.adaptive = true;
      V.scale = 55; V.rotFrame = true;
      const sun = makeBody(1, 0, 0, 0, 0, "#ffd166", "Sun");
      bodies.push(sun);
      const aJ = 5.2026, mJ = 9.546e-4;
      const vJ = Math.sqrt(Gv() * (1 + mJ) / aJ);
      bodies.push(makeBody(mJ, aJ, 0, 0, vJ, "#e0b080", "Jupiter"));
      zeroTotalMomentum();
      // Trojan swarms ±60° with small offsets, co-rotating circular speed
      const rnd = mulberry32(rngSeed + 1);
      for (const sgn of [1, -1]) {
        for (let i = 0; i < 5; i++) {
          const ang = sgn * Math.PI / 3 + (rnd() - 0.5) * 0.12;
          const r = aJ * (1 + (rnd() - 0.5) * 0.02);
          const v = Math.sqrt(Gv() * 1 / r);
          bodies.push(makeBody(1e-12,
            r * Math.cos(ang), r * Math.sin(ang),
            -v * Math.sin(ang), v * Math.cos(ang),
            sgn > 0 ? "#5ce8c5" : "#c792ea",
            (sgn > 0 ? "L4" : "L5") + ` trojan ${i + 1}`));
        }
      }
    },
  },

  binary: {
    label: "Binary stars + circumbinary planet",
    info: "Two suns orbiting their barycenter with a planet around the pair.",
    build() {
      P.dt = 0.2; P.sub = 12; P.eps = 0;
      P.integrator = "yoshida4"; P.adaptive = true;
      V.scale = 250;
      const d = 0.5;
      const v = Math.sqrt(Gv() * 1 / (2 * d));
      bodies.push(makeBody(1, -d / 2, 0, 0, -v, "#ffd166", "star A"));
      bodies.push(makeBody(1,  d / 2, 0, 0,  v, "#ff6b8a", "star B"));
      const rp = 1.6, vp = Math.sqrt(Gv() * 2 / rp);
      bodies.push(makeBody(3e-6, 0, rp, vp, 0, "#7fd4ff", "planet"));
    },
  },

  figure8: {
    label: "Three-body figure-8 (1 M☉ each)",
    info: "The Chenciner–Montgomery choreography with solar-mass stars at AU scale; " +
          "exact initial conditions, period ≈ 368 days.",
    build() {
      P.dt = 0.1; P.sub = 30; P.eps = 0;
      P.integrator = "yoshida4"; P.adaptive = true;
      V.scale = 280;
      const vs = Math.sqrt(Gv() * 1 / 1);   // velocity scale for G=G☉, m=1, L=1 AU
      const p = [[-0.97000436, 0.24308753], [0.97000436, -0.24308753], [0, 0]];
      const v3 = [-0.93240737, -0.86473146];
      const vel = [[-v3[0] / 2, -v3[1] / 2], [-v3[0] / 2, -v3[1] / 2], v3];
      for (let i = 0; i < 3; i++) {
        bodies.push(makeBody(1, p[i][0], p[i][1], vel[i][0] * vs, vel[i][1] * vs,
          undefined, `star ${i + 1}`));
      }
    },
  },

  accretion: {
    label: "Planetesimal accretion (seeded chaos)",
    info: "A young star with 16 planetesimals on crossing orbits, merging on collision. " +
          "Seeded RNG: the same seed always reproduces the same system.",
    build() {
      P.dt = 0.5; P.sub = 12; P.eps = 0.005;
      P.integrator = "verlet"; P.adaptive = true; P.merge = true;
      V.scale = 120;
      const sun = makeBody(1, 0, 0, 0, 0, "#ffd166", "star");
      bodies.push(sun);
      const rnd = mulberry32(rngSeed);
      for (let i = 0; i < 16; i++) {
        const a = rnd() * 2 * Math.PI;
        const r = 0.6 + rnd() * 2.6;
        const m = 1e-7 * Math.pow(10, rnd() * 2);     // 1e-7 .. 1e-5 M☉
        const v = Math.sqrt(Gv() / r) * (0.85 + rnd() * 0.3);
        bodies.push(makeBody(m, r * Math.cos(a), r * Math.sin(a),
          -v * Math.sin(a), v * Math.cos(a), undefined, `pl-${i + 1}`));
      }
      zeroTotalMomentum();
    },
  },

  empty: {
    label: "Empty space (build your own)",
    info: "Use 🚀 Launch mode to place bodies, or edit them numerically in the panel.",
    build() {
      P.dt = 0.5; P.sub = 12; P.eps = 0;
      P.integrator = "verlet"; P.adaptive = true;
      V.scale = 120;
    },
  },
};

function zeroTotalMomentum() {
  let M = 0, px = 0, py = 0;
  for (const b of bodies) { M += b.m; px += b.m * b.vx; py += b.m * b.vy; }
  if (M > 0) for (const b of bodies) { b.vx -= px / M; b.vy -= py / M; }
}
