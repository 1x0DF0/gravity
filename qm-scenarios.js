/* ============================================================
   Quantum Lab — scenarios
   Internal units ħ = m = 1 (each scenario states its real-unit
   mapping). Every scenario carries its own analytic prediction
   and measures the simulation against it, live.
   ============================================================ */

"use strict";

// Airy-function zeros: bouncer energies E_n = (ħ²mg²/2)^⅓ · aₙ
const AIRY_ZEROS = [2.33810741, 4.08794944, 5.52055983, 6.78670809];
// real neutron numbers for the qBOUNCE mapping
const BOUNCER_E1_PEV = 1.4074;            // (ħ²mg²/2)^⅓·a₁ for a neutron, in peV
const BOUNCER_LG_UM = 5.87;               // (ħ²/2m²g)^⅓ in µm

const QM_SCENARIOS = {

  free: {
    label: "Free wave packet (dispersion)",
    info: "A free Gaussian packet spreads as σ(t) = σ₀√(1+(ħt/2mσ₀²)²) — pure matter-wave " +
          "dispersion. The dashed envelope is the analytic prediction.",
    dim: 1, n: 1024, L: 400, dt: 0.05, perFrame: 4,
    params: [
      { key: "sigma", label: "initial width σ₀", min: 2, max: 12, step: 0.1, val: 5 },
      { key: "k0", label: "momentum k₀", min: 0, max: 2, step: 0.05, val: 0.8 },
    ],
    build(S, p) {
      S.setPotential(() => 0);
      S.gaussian({ x0: -100, sx: p.sigma, kx: p.k0 });
      S.normalize();
      this._hist = [];
    },
    sigmaPred(p, t) { return p.sigma * Math.sqrt(1 + (t / (2 * p.sigma * p.sigma)) ** 2); },
    tick(S, p) {
      const m = S.moments();
      this._hist.push([S.t, Math.sqrt(m.varX), this.sigmaPred(p, S.t)]);
      if (this._hist.length > 500) this._hist.shift();
    },
    measure(S, p) {
      const m = S.moments();
      const sm = Math.sqrt(m.varX), sp = this.sigmaPred(p, S.t);
      return [
        `σₓ measured  = ${sm.toFixed(3)}`,
        `σₓ analytic  = ${sp.toFixed(3)}   (${(100 * (sm - sp) / sp).toFixed(2)}% off)`,
        `⟨x⟩ = ${m.x.toFixed(2)}  (classical: ${(-100 + p.k0 * S.t).toFixed(2)})`,
      ];
    },
    graph: { title: "σₓ(t): measured vs analytic", series: 2 },
  },

  harmonic: {
    label: "Harmonic oscillator (coherent state)",
    info: "A displaced ground state (coherent state) oscillates at the classical period " +
          "2π/ω without spreading — the closest quantum mechanics gets to classical motion. " +
          "Eigenenergies below are found by imaginary-time relaxation; exact: Eₙ = (n+½)ω.",
    dim: 1, n: 1024, L: 120, dt: 0.01, perFrame: 8,
    params: [
      { key: "omega", label: "trap frequency ω", min: 0.2, max: 1.5, step: 0.01, val: 0.5 },
      { key: "x0", label: "displacement x₀", min: 0, max: 25, step: 0.5, val: 12 },
    ],
    build(S, p) {
      S.setPotential(x => 0.5 * p.omega * p.omega * x * x);
      S.gaussian({ x0: p.x0, sx: Math.sqrt(0.5 / p.omega) });
      S.normalize();
      this._hist = [];
    },
    eigen: { count: 4, exact: (n, p) => (n + 0.5) * p.omega, unit: "ħω₀" },
    tick(S, p) {
      const m = S.moments();
      this._hist.push([S.t, m.x, p.x0 * Math.cos(p.omega * S.t)]);
      if (this._hist.length > 500) this._hist.shift();
    },
    measure(S, p) {
      const m = S.moments();
      return [
        `⟨x⟩ measured  = ${m.x.toFixed(3)}`,
        `⟨x⟩ classical = ${(p.x0 * Math.cos(p.omega * S.t)).toFixed(3)}`,
        `σₓ = ${Math.sqrt(m.varX).toFixed(3)}  (coherent: ${Math.sqrt(0.5 / p.omega).toFixed(3)}, constant)`,
        `period 2π/ω = ${(2 * Math.PI / p.omega).toFixed(2)}`,
      ];
    },
    graph: { title: "⟨x⟩(t): quantum vs classical", series: 2 },
  },

  tunnel: {
    label: "Tunneling through a barrier",
    info: "The packet's mean energy is below the barrier top, yet part of it gets through. " +
          "Measured transmission is compared with the analytic plane-wave result averaged " +
          "over the packet's momentum distribution: T = ∫|φ(k)|²T(k)dk.",
    dim: 1, n: 2048, L: 600, dt: 0.04, perFrame: 6,
    params: [
      { key: "V0", label: "barrier height V₀", min: 0.2, max: 2, step: 0.01, val: 1.0 },
      { key: "a", label: "barrier width a", min: 0.5, max: 6, step: 0.1, val: 2 },
      { key: "k0", label: "packet momentum k₀", min: 0.5, max: 2, step: 0.01, val: 1.2 },
    ],
    build(S, p) {
      S.setPotential(x => (x > -p.a / 2 && x < p.a / 2) ? p.V0 : 0);
      S.gaussian({ x0: -90, sx: 10, kx: p.k0 });
      S.normalize();
      this._hist = [];
    },
    // the prediction uses the barrier as actually discretized on the grid
    // (T depends exponentially on width, so the nominal width is not honest)
    effWidth(S, p) {
      let c = 0;
      for (let i = 0; i < S.n; i++) if (S.V[i] > p.V0 / 2) c++;
      return c * S.dx;
    },
    tick(S, p) {
      this._hist.push([S.t, S.probIn(p.a / 2 + 5, 1e9), packetT(p.k0, 10, p.V0, this.effWidth(S, p))]);
      if (this._hist.length > 500) this._hist.shift();
    },
    measure(S, p) {
      const T = S.probIn(p.a / 2 + 5, 1e9);
      const Tpred = packetT(p.k0, 10, p.V0, this.effWidth(S, p));
      const E = p.k0 * p.k0 / 2;
      return [
        `⟨E⟩ = ${E.toFixed(3)}   V₀ = ${p.V0.toFixed(2)}  (E/V₀ = ${(E / p.V0).toFixed(2)})`,
        `T measured  = ${T.toFixed(4)}  (settles when packet clears)`,
        `T analytic  = ${Tpred.toFixed(4)}`,
      ];
    },
    graph: { title: "transmitted probability vs analytic T", series: 2 },
  },

  bouncer: {
    label: "Quantum bouncer (qBOUNCE, real neutrons)",
    info: "A neutron bouncing on a mirror in Earth's gravity occupies discrete quantum " +
          "states — measured by Nesvizhevsky et al. (2002) and qBOUNCE. Eigenvalues below " +
          "are Airy-function zeros; the real-unit column uses the neutron mass: " +
          "E₁ = 1.41 peV, length scale 5.87 µm. The wall at x = 0 is exact " +
          "(odd-parity sector of V = g|x|).",
    dim: 1, n: 1024, L: 60, dt: 0.002, perFrame: 20,
    halfDomain: true,        // x<0 is the odd-parity mirror image, not physics
    params: [
      { key: "h0", label: "drop height", min: 2, max: 12, step: 0.1, val: 6 },
    ],
    build(S, p) {
      S.setPotential(x => Math.abs(x));        // g = 1; wall via odd projection
      S.oddProject = true;
      S.gaussian({ x0: p.h0, sx: 1.2 });
      S.applyOdd();
      S.normalize();
      this._hist = [];
    },
    eigen: {
      count: 4,
      exact: n => Math.cbrt(0.5) * AIRY_ZEROS[n],
      unit: "ε₀",
      real: E => `${(E / (Math.cbrt(0.5) * AIRY_ZEROS[0]) * BOUNCER_E1_PEV).toFixed(3)} peV`,
    },
    tick(S) {
      const m = S.moments(0);
      this._hist.push([S.t, m.x, null]);
      if (this._hist.length > 500) this._hist.shift();
    },
    measure(S, p) {
      const m = S.moments(0);
      return [
        `⟨x⟩ = ${m.x.toFixed(2)}  (height in units of ℓ_g = ${BOUNCER_LG_UM} µm for a neutron)`,
        `classical bounce period 2√(2h₀/g) = ${(2 * Math.sqrt(2 * p.h0)).toFixed(2)}`,
        `collapses & revivals: the packet dephases, then partially reassembles`,
      ];
    },
    graph: { title: "⟨x⟩(t): bounce, collapse, revival", series: 1 },
  },

  gravimeter: {
    label: "Atom-interferometer gravimeter",
    info: "A matter wave held in a superposition of two heights accumulates relative phase " +
          "Δφ(t) = mgΔh·t/ħ — the gravitationally induced quantum phase first seen in the " +
          "COW neutron experiment (1975); Kasevich–Chu atom interferometers measure g this " +
          "way to ~10⁻⁹. Watch the measured phase slope match mgΔh exactly.",
    dim: 1, n: 1024, L: 120, dt: 0.01, perFrame: 8,
    params: [
      { key: "g", label: "gravity g", min: 0, max: 0.12, step: 0.002, val: 0.05 },
      { key: "dh", label: "height separation Δh", min: 8, max: 40, step: 1, val: 20 },
    ],
    build(S, p) {
      const w = 1;                              // trap frequency
      const h1 = -p.dh / 2, h2 = p.dh / 2;
      S.setPotential(x =>
        0.5 * w * w * Math.min((x - h1) ** 2, (x - h2) ** 2) + p.g * x);
      const s0 = Math.sqrt(0.5 / w);
      const xeq = p.g / (w * w);                // tilt shifts both minima equally
      S.gaussian({ x0: h1 - xeq, sx: s0, amp: 1 });
      S.gaussian({ x0: h2 - xeq, sx: s0, amp: 1, add: true });
      S.normalize();
      this._h1 = h1 - xeq; this._h2 = h2 - xeq;
      this._prev = null; this._acc = 0; this._hist = [];
    },
    phaseAt(S, x0) {
      const i = Math.round((x0 + S.L / 2) / S.dx);
      return Math.atan2(S.im[i], S.re[i]);
    },
    tick(S, p) {
      let d = this.phaseAt(S, this._h1) - this.phaseAt(S, this._h2);
      if (this._prev !== null) {
        let dd = d - this._prev;
        while (dd > Math.PI) dd -= 2 * Math.PI;
        while (dd < -Math.PI) dd += 2 * Math.PI;
        this._acc += dd;
      }
      this._prev = d;
      this._hist.push([S.t, this._acc, p.g * p.dh * S.t]);
      if (this._hist.length > 500) this._hist.shift();
    },
    measure(S, p) {
      const slope = S.t > 1 ? this._acc / S.t : 0;
      return [
        `Δφ accumulated = ${this._acc.toFixed(3)} rad`,
        `dΔφ/dt measured  = ${slope.toFixed(5)} rad/t`,
        `dΔφ/dt = mgΔh/ħ  = ${(p.g * p.dh).toFixed(5)} rad/t  (predicted)`,
        slope !== 0 ? `g inferred from phase = ${(slope / p.dh).toFixed(5)}  (true: ${p.g.toFixed(5)})` : "",
      ];
    },
    graph: { title: "Δφ(t): measured vs mgΔh·t/ħ", series: 2 },
  },

  doubleslit: {
    label: "Double slit (matter waves)",
    info: "A packet with de Broglie wavelength λ = 2π/k₀ passes two slits; the arrival " +
          "pattern accumulates on the screen (right strip). Fringe spacing is compared " +
          "with λD/d. Edges absorb (norm intentionally decays once waves leave).",
    dim: 2, nx: 256, ny: 256, Lx: 120, Ly: 120, dt: 0.02, perFrame: 3,
    params: [
      { key: "k0", label: "momentum k₀", min: 1, max: 3, step: 0.05, val: 2 },
      { key: "d", label: "slit separation d", min: 6, max: 20, step: 0.5, val: 10 },
      { key: "sw", label: "slit width", min: 1.5, max: 6, step: 0.25, val: 3 },
    ],
    screenX: 40,
    build(S, p) {
      const wallW = 1.5, V0 = 60;
      S.setPotential((x, y) => {
        if (Math.abs(x) > wallW / 2) return 0;
        const inSlit = Math.abs(Math.abs(y) - p.d / 2) < p.sw / 2;
        return inSlit ? 0 : V0;
      });
      S.setAbsorber(0.08, 0.35);
      S.gaussian({ x0: -32, y0: 0, sx: 7, sy: 18, kx: p.k0 });
      S.normalize();
      this._screen = new Float64Array(S.ny);
    },
    tick(S) {
      const ix = Math.round((this.screenX + S.Lx / 2) / S.dx);
      for (let j = 0; j < S.ny; j++) {
        const idx = j * S.nx + ix;
        this._screen[j] += S.re[idx] ** 2 + S.im[idx] ** 2;
      }
    },
    fringeSpacing(S) {
      // textbook estimator: distance between the first minima either side of
      // the central maximum — equals λD/d exactly, even at finite D where the
      // off-axis spacing chirps (path difference is d·sinθ but y = D·tanθ)
      const n = S.ny;
      const s = Float64Array.from(this._screen);
      for (let pass = 0; pass < 2; pass++) {     // light smoothing
        const c = s.slice();
        for (let i = 1; i < n - 1; i++) s[i] = (c[i - 1] + 2 * c[i] + c[i + 1]) / 4;
      }
      const c0 = n / 2;
      const firstMin = dir => {
        for (let i = c0 + 2 * dir; i > 2 && i < n - 2; i += dir) {
          if (s[i] < s[i - 1] && s[i] <= s[i + 1]) {
            // parabolic sub-cell refinement of the minimum position
            const denom = s[i - 1] - 2 * s[i] + s[i + 1];
            const off = Math.abs(denom) > 1e-30
              ? Math.max(-0.5, Math.min(0.5, 0.5 * (s[i - 1] - s[i + 1]) / denom)) : 0;
            return i + off;
          }
        }
        return null;
      };
      const lo = firstMin(-1), hi = firstMin(1);
      return (lo !== null && hi !== null) ? (hi - lo) * S.dy : null;
    },
    measure(S, p) {
      const lam = 2 * Math.PI / p.k0;
      const pred = lam * this.screenX / p.d;
      const tot = this._screen.reduce((a, b) => a + b, 0);
      const meas = tot > 1 ? this.fringeSpacing(S) : null;
      return [
        `λ = 2π/k₀ = ${lam.toFixed(3)}   D = ${this.screenX}   d = ${p.d}`,
        `fringe spacing predicted λD/d = ${pred.toFixed(2)}`,
        meas ? `fringe spacing measured = ${meas.toFixed(2)}  (${(100 * (meas - pred) / pred).toFixed(1)}% off)`
             : `fringe spacing measured = … (accumulating)`,
      ];
    },
    graph: { title: "screen arrival pattern", screen: true },
  },

  orbit: {
    label: "Ehrenfest: packet orbiting a 1/r well",
    info: "A wave packet in an attractive Coulomb-like well, with the classical trajectory " +
          "(white) integrated from the same initial ⟨x⟩, ⟨p⟩. Ehrenfest's theorem: ⟨x⟩ " +
          "follows the classical path — until the packet spreads around the nucleus and " +
          "the correspondence visibly fails. This is Gravity Lab's orbit, quantized.",
    dim: 2, nx: 256, ny: 256, Lx: 140, Ly: 140, dt: 0.02, perFrame: 4,
    params: [
      { key: "alpha", label: "well strength α", min: 1, max: 10, step: 0.1, val: 5 },
      { key: "r0", label: "initial radius r₀", min: 12, max: 40, step: 1, val: 22 },
    ],
    build(S, p) {
      const s = 2.5;
      S.setPotential((x, y) => -p.alpha / Math.sqrt(x * x + y * y + s * s));
      S.setAbsorber(0.06, 0.25);
      const r0 = p.r0;
      const aMag = p.alpha * r0 / Math.pow(r0 * r0 + s * s, 1.5);
      const v = Math.sqrt(r0 * aMag);          // circular-ish speed
      S.gaussian({ x0: r0, y0: 0, sx: 3.5, sy: 3.5, kx: 0, ky: v });
      S.normalize();
      this._cl = { x: r0, y: 0, vx: 0, vy: v, s, trail: [] };
      this._hist = [];
    },
    clAccel(x, y, p) {
      const s = this._cl.s;
      const f = -p.alpha / Math.pow(x * x + y * y + s * s, 1.5);
      return [f * x, f * y];
    },
    tick(S, p, dtFrame) {
      const c = this._cl;                       // RK4, same elapsed time
      const h = dtFrame;
      const d = (x, y, vx, vy) => { const [ax, ay] = this.clAccel(x, y, p); return [vx, vy, ax, ay]; };
      const k1 = d(c.x, c.y, c.vx, c.vy);
      const k2 = d(c.x + h / 2 * k1[0], c.y + h / 2 * k1[1], c.vx + h / 2 * k1[2], c.vy + h / 2 * k1[3]);
      const k3 = d(c.x + h / 2 * k2[0], c.y + h / 2 * k2[1], c.vx + h / 2 * k2[2], c.vy + h / 2 * k2[3]);
      const k4 = d(c.x + h * k3[0], c.y + h * k3[1], c.vx + h * k3[2], c.vy + h * k3[3]);
      c.x += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
      c.y += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
      c.vx += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
      c.vy += h / 6 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
      c.trail.push(c.x, c.y);
      if (c.trail.length > 700) c.trail.splice(0, c.trail.length - 700);
      const m = S.moments();
      this._hist.push([S.t, Math.hypot(m.x - c.x, m.y - c.y), null]);
      if (this._hist.length > 500) this._hist.shift();
    },
    measure(S) {
      const m = S.moments();
      const c = this._cl;
      return [
        `⟨r⟩ quantum  = (${m.x.toFixed(1)}, ${m.y.toFixed(1)})`,
        `r classical = (${c.x.toFixed(1)}, ${c.y.toFixed(1)})`,
        `|⟨r⟩ − r_cl| = ${Math.hypot(m.x - c.x, m.y - c.y).toFixed(2)}   σₓ = ${Math.sqrt(m.varX).toFixed(2)}`,
      ];
    },
    graph: { title: "|⟨r⟩ − r_classical|(t): Ehrenfest breakdown", series: 1 },
  },

  entangle: {
    label: "Two particles: entanglement from collision",
    info: "This plot is NOT physical space — it is the configuration space ψ(x₁,x₂) of TWO " +
          "particles on a line. They start as a product state (purity 1), scatter off each " +
          "other, and become entangled: the state no longer factorizes and the Schmidt " +
          "number K = 1/Tr ρ₁² grows. This exponential growth of state space with particle " +
          "count is Feynman's argument for quantum computers.",
    dim: 2, nx: 128, ny: 128, Lx: 100, Ly: 100, dt: 0.02, perFrame: 4,
    params: [
      // at E_rel = k² = V₀ the collision splits 50/50 → maximal K ≈ 2
      { key: "V0", label: "interaction strength V₀", min: 0, max: 6, step: 0.1, val: 2.5 },
      { key: "k", label: "particle momentum", min: 0.3, max: 2.5, step: 0.02, val: 1.58 },
    ],
    axes: ["x₁ (particle 1)", "x₂ (particle 2)"],
    build(S, p) {
      const w = 1.8;
      S.setPotential((x1, x2) => p.V0 * Math.exp(-((x1 - x2) ** 2) / (2 * w * w)));
      S.setAbsorber(0.07, 0.3);
      S.gaussian({ x0: -20, y0: 20, sx: 5, sy: 5, kx: p.k, ky: -p.k });
      S.normalize();
      this._hist = [];
      this._K = 1; this._ctr = 0;
    },
    tick(S) {
      if (this._ctr++ % 6 === 0) this._K = 1 / S.purity();
      this._hist.push([S.t, this._K, null]);
      if (this._hist.length > 500) this._hist.shift();
    },
    measure(S, p) {
      return [
        `Schmidt number K = 1/Tr ρ₁² = ${this._K.toFixed(3)}`,
        `K = 1 ⟺ product state (no entanglement)`,
        `K → 2 at a 50/50 collision (E_rel = k² = V₀, i.e. k = ${Math.sqrt(p.V0).toFixed(2)}):`,
        `the two outcomes (bounced / passed through) are the Schmidt modes`,
      ];
    },
    graph: { title: "Schmidt number K(t)", series: 1 },
  },
};

if (typeof module !== "undefined") {
  module.exports = { QM_SCENARIOS, AIRY_ZEROS, BOUNCER_E1_PEV };
}
