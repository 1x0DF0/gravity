/* ============================================================
   Quantum Lab — TDSE engine
   Solves  iħ ∂ψ/∂t = [ −ħ²/2m ∇² + V(x[,y]) ] ψ   (ħ = 1 internal)
   with the split-step Fourier method: exactly unitary, spectral
   accuracy. Pure simulation logic, no DOM.
   ============================================================ */

"use strict";

// ---------------- FFT (iterative radix-2, separate re/im) ----------------

const _fftTables = new Map();
function fftTables(n) {
  let t = _fftTables.get(n);
  if (t) return t;
  const bits = Math.round(Math.log2(n));
  const rev = new Uint32Array(n);
  for (let i = 0; i < n; i++) {
    let r = 0, x = i;
    for (let b = 0; b < bits; b++) { r = (r << 1) | (x & 1); x >>= 1; }
    rev[i] = r;
  }
  const cos = new Float64Array(n / 2), sin = new Float64Array(n / 2);
  for (let i = 0; i < n / 2; i++) {
    cos[i] = Math.cos(-2 * Math.PI * i / n);
    sin[i] = Math.sin(-2 * Math.PI * i / n);
  }
  t = { rev, cos, sin };
  _fftTables.set(n, t);
  return t;
}

function fft(re, im, inv) {
  const n = re.length;
  const { rev, cos, sin } = fftTables(n);
  for (let i = 0; i < n; i++) {
    const j = rev[i];
    if (j > i) {
      let t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const half = len >> 1, step = n / len;
    for (let i = 0; i < n; i += len) {
      for (let j = 0; j < half; j++) {
        const k = j * step;
        const wr = cos[k], wi = inv ? -sin[k] : sin[k];
        const a = i + j, b = a + half;
        const tr = re[b] * wr - im[b] * wi;
        const ti = re[b] * wi + im[b] * wr;
        re[b] = re[a] - tr; im[b] = im[a] - ti;
        re[a] += tr; im[a] += ti;
      }
    }
  }
  if (inv) {
    const s = 1 / n;
    for (let i = 0; i < n; i++) { re[i] *= s; im[i] *= s; }
  }
}

// spatial frequency k for index j on an n-point grid of length L
function kOf(j, n, L) {
  return 2 * Math.PI * (j < n / 2 ? j : j - n) / L;
}

// ---------------- quantum state (1D or 2D) ----------------

class QMState {
  constructor(opts) {
    this.dim = opts.dim || 1;
    this.m = opts.m || 1;
    if (this.dim === 1) {
      this.n = opts.n;
      this.L = opts.L;
      this.dx = this.L / this.n;
      this.size = this.n;
    } else {
      this.nx = opts.nx; this.ny = opts.ny;
      this.Lx = opts.Lx; this.Ly = opts.Ly;
      this.dx = this.Lx / this.nx; this.dy = this.Ly / this.ny;
      this.size = this.nx * this.ny;
    }
    this.re = new Float64Array(this.size);
    this.im = new Float64Array(this.size);
    this.V = new Float64Array(this.size);
    this.mask = null;           // absorbing-boundary factors (optional)
    this.oddProject = false;    // restrict to odd sector (hard wall at x = 0)
    this.t = 0;
    this._prep = null;          // cached propagator phases
    this._sr = new Float64Array(this.dim === 2 ? Math.max(this.nx, this.ny) : 0);
    this._si = new Float64Array(this.dim === 2 ? Math.max(this.nx, this.ny) : 0);
  }

  x(i) { return -((this.dim === 1 ? this.L : this.Lx) / 2) + i * this.dx; }
  y(j) { return -(this.Ly / 2) + j * this.dy; }
  cellMeasure() { return this.dim === 1 ? this.dx : this.dx * this.dy; }

  setPotential(fn) {
    if (this.dim === 1) {
      for (let i = 0; i < this.n; i++) this.V[i] = fn(this.x(i));
    } else {
      for (let j = 0; j < this.ny; j++)
        for (let i = 0; i < this.nx; i++)
          this.V[j * this.nx + i] = fn(this.x(i), this.y(j));
    }
    this._prep = null;
  }

  // Gaussian packet  exp(−(r−r0)²/4σ² + i k·r); multiplies onto existing ψ if add=true
  gaussian(opts) {
    const { x0 = 0, y0 = 0, sx = 5, sy = 5, kx = 0, ky = 0, amp = 1, add = false } = opts;
    if (this.dim === 1) {
      for (let i = 0; i < this.n; i++) {
        const x = this.x(i) - x0;
        const g = amp * Math.exp(-x * x / (4 * sx * sx));
        const ph = kx * this.x(i);
        const r = g * Math.cos(ph), im_ = g * Math.sin(ph);
        if (add) { this.re[i] += r; this.im[i] += im_; }
        else { this.re[i] = r; this.im[i] = im_; }
      }
    } else {
      for (let j = 0; j < this.ny; j++) {
        for (let i = 0; i < this.nx; i++) {
          const x = this.x(i) - x0, y = this.y(j) - y0;
          const g = amp * Math.exp(-x * x / (4 * sx * sx) - y * y / (4 * sy * sy));
          const ph = kx * this.x(i) + ky * this.y(j);
          const idx = j * this.nx + i;
          const r = g * Math.cos(ph), im_ = g * Math.sin(ph);
          if (add) { this.re[idx] += r; this.im[idx] += im_; }
          else { this.re[idx] = r; this.im[idx] = im_; }
        }
      }
    }
  }

  norm() {
    let s = 0;
    for (let i = 0; i < this.size; i++) s += this.re[i] * this.re[i] + this.im[i] * this.im[i];
    return s * this.cellMeasure();
  }

  normalize() {
    const s = 1 / Math.sqrt(this.norm());
    for (let i = 0; i < this.size; i++) { this.re[i] *= s; this.im[i] *= s; }
  }

  applyOdd() {     // ψ(x) ← [ψ(x) − ψ(−x)]/2 — exact hard wall at x = 0 (1D)
    const n = this.n;
    for (let i = 1; i < n / 2; i++) {
      const j = n - i;
      const r = 0.5 * (this.re[i] - this.re[j]), m = 0.5 * (this.im[i] - this.im[j]);
      this.re[i] = r; this.im[i] = m;
      this.re[j] = -r; this.im[j] = -m;
    }
    this.re[0] = this.im[0] = 0;
    this.re[n / 2] = this.im[n / 2] = 0;
  }

  // precompute exp(−iV dt/2) and exp(−i k²/2m dt) for the current dt
  prepare(dt) {
    if (this._prep && this._prep.dt === dt) return this._prep;
    const p = { dt };
    p.vc = new Float64Array(this.size); p.vs = new Float64Array(this.size);
    for (let i = 0; i < this.size; i++) {
      const ph = -this.V[i] * dt / 2;
      p.vc[i] = Math.cos(ph); p.vs[i] = Math.sin(ph);
    }
    p.kc = new Float64Array(this.size); p.ks = new Float64Array(this.size);
    if (this.dim === 1) {
      for (let i = 0; i < this.n; i++) {
        const k = kOf(i, this.n, this.L);
        const ph = -k * k / (2 * this.m) * dt;
        p.kc[i] = Math.cos(ph); p.ks[i] = Math.sin(ph);
      }
    } else {
      for (let j = 0; j < this.ny; j++) {
        const ky = kOf(j, this.ny, this.Ly);
        for (let i = 0; i < this.nx; i++) {
          const kx = kOf(i, this.nx, this.Lx);
          const ph = -(kx * kx + ky * ky) / (2 * this.m) * dt;
          const idx = j * this.nx + i;
          p.kc[idx] = Math.cos(ph); p.ks[idx] = Math.sin(ph);
        }
      }
    }
    this._prep = p;
    return p;
  }

  _mulPhase(c, s) {
    for (let i = 0; i < this.size; i++) {
      const r = this.re[i], m = this.im[i];
      this.re[i] = r * c[i] - m * s[i];
      this.im[i] = r * s[i] + m * c[i];
    }
  }

  fftAll(inv) {
    if (this.dim === 1) { fft(this.re, this.im, inv); return; }
    const { nx, ny } = this;
    for (let j = 0; j < ny; j++) {
      fft(this.re.subarray(j * nx, j * nx + nx), this.im.subarray(j * nx, j * nx + nx), inv);
    }
    const sr = this._sr, si = this._si;
    for (let i = 0; i < nx; i++) {
      for (let j = 0; j < ny; j++) { sr[j] = this.re[j * nx + i]; si[j] = this.im[j * nx + i]; }
      fft(sr.subarray(0, ny), si.subarray(0, ny), inv);
      for (let j = 0; j < ny; j++) { this.re[j * nx + i] = sr[j]; this.im[j * nx + i] = si[j]; }
    }
  }

  // one split-step:  e^{−iVdt/2} · F⁻¹ e^{−ik²dt/2m} F · e^{−iVdt/2}
  step(dt, steps = 1) {
    const p = this.prepare(dt);
    for (let s = 0; s < steps; s++) {
      this._mulPhase(p.vc, p.vs);
      this.fftAll(false);
      this._mulPhase(p.kc, p.ks);
      this.fftAll(true);
      this._mulPhase(p.vc, p.vs);
      if (this.mask) {
        for (let i = 0; i < this.size; i++) { this.re[i] *= this.mask[i]; this.im[i] *= this.mask[i]; }
      }
      if (this.oddProject) this.applyOdd();
      this.t += dt;
    }
  }

  // imaginary-time step (for eigenstates): e^{−V dτ/2} F⁻¹ e^{−k²dτ/2m} F e^{−V dτ/2}
  imagStep(dtau) {
    if (!this._iprep || this._iprep.dtau !== dtau) {
      const ip = { dtau, v: new Float64Array(this.size), k: new Float64Array(this.size) };
      for (let i = 0; i < this.size; i++) ip.v[i] = Math.exp(-this.V[i] * dtau / 2);
      if (this.dim === 1) {
        for (let i = 0; i < this.n; i++) {
          const k = kOf(i, this.n, this.L);
          ip.k[i] = Math.exp(-k * k / (2 * this.m) * dtau);
        }
      } else {
        for (let j = 0; j < this.ny; j++) {
          const ky = kOf(j, this.ny, this.Ly);
          for (let i = 0; i < this.nx; i++) {
            const kx = kOf(i, this.nx, this.Lx);
            ip.k[j * this.nx + i] = Math.exp(-(kx * kx + ky * ky) / (2 * this.m) * dtau);
          }
        }
      }
      this._iprep = ip;
    }
    const ip = this._iprep;
    for (let i = 0; i < this.size; i++) { this.re[i] *= ip.v[i]; this.im[i] *= ip.v[i]; }
    this.fftAll(false);
    for (let i = 0; i < this.size; i++) { this.re[i] *= ip.k[i]; this.im[i] *= ip.k[i]; }
    this.fftAll(true);
    for (let i = 0; i < this.size; i++) { this.re[i] *= ip.v[i]; this.im[i] *= ip.v[i]; }
  }

  // ⟨φ|ψ⟩ with grid measure
  dot(other) {
    let r = 0, m = 0;
    for (let i = 0; i < this.size; i++) {
      r += other.re[i] * this.re[i] + other.im[i] * this.im[i];
      m += other.re[i] * this.im[i] - other.im[i] * this.re[i];
    }
    const c = this.cellMeasure();
    return [r * c, m * c];
  }

  subtractProjection(other) {   // ψ ← ψ − ⟨φ|ψ⟩φ
    const [r, m] = this.dot(other);
    for (let i = 0; i < this.size; i++) {
      this.re[i] -= r * other.re[i] - m * other.im[i];
      this.im[i] -= r * other.im[i] + m * other.re[i];
    }
  }

  // ---------------- observables ----------------

  // ⟨x⟩, ⟨x²⟩ (and y in 2D) over |ψ|², normalized. xMin restricts to the
  // physical region when a mirror-trick half-domain is in use (1D only).
  moments(xMin = -Infinity) {
    let W = 0, mx = 0, mx2 = 0, my = 0, my2 = 0;
    if (this.dim === 1) {
      for (let i = 0; i < this.n; i++) {
        const x = this.x(i);
        if (x < xMin) continue;
        const w = this.re[i] * this.re[i] + this.im[i] * this.im[i];
        W += w; mx += w * x; mx2 += w * x * x;
      }
    } else {
      for (let j = 0; j < this.ny; j++) {
        const y = this.y(j);
        for (let i = 0; i < this.nx; i++) {
          const w = this.re[j * this.nx + i] ** 2 + this.im[j * this.nx + i] ** 2;
          const x = this.x(i);
          W += w; mx += w * x; mx2 += w * x * x; my += w * y; my2 += w * y * y;
        }
      }
    }
    mx /= W; mx2 /= W; my /= W; my2 /= W;
    return { x: mx, y: my, varX: Math.max(0, mx2 - mx * mx), varY: Math.max(0, my2 - my * my) };
  }

  // momentum moments + kinetic energy via the spectrum (needs a scratch copy)
  spectral() {
    const re = this.re.slice(), im = this.im.slice();
    const save = [this.re, this.im];
    this.re = re; this.im = im;
    this.fftAll(false);
    let W = 0, mk = 0, mk2 = 0, mky = 0, mky2 = 0, T = 0;
    if (this.dim === 1) {
      for (let i = 0; i < this.n; i++) {
        const w = re[i] * re[i] + im[i] * im[i];
        const k = kOf(i, this.n, this.L);
        W += w; mk += w * k; mk2 += w * k * k;
      }
      T = mk2 / (2 * this.m);
    } else {
      for (let j = 0; j < this.ny; j++) {
        const ky = kOf(j, this.ny, this.Ly);
        for (let i = 0; i < this.nx; i++) {
          const kx = kOf(i, this.nx, this.Lx);
          const w = re[j * this.nx + i] ** 2 + im[j * this.nx + i] ** 2;
          W += w; mk += w * kx; mk2 += w * kx * kx; mky += w * ky; mky2 += w * ky * ky;
        }
      }
      T = (mk2 + mky2) / (2 * this.m);
    }
    [this.re, this.im] = save;
    mk /= W; mk2 /= W; mky /= W; mky2 /= W; T /= W;
    return { px: mk, py: mky, varPx: Math.max(0, mk2 - mk * mk),
             varPy: Math.max(0, mky2 - mky * mky), T };
  }

  potentialEnergy() {
    let W = 0, EV = 0;
    for (let i = 0; i < this.size; i++) {
      const w = this.re[i] * this.re[i] + this.im[i] * this.im[i];
      W += w; EV += w * this.V[i];
    }
    return EV / W;
  }

  energy() {
    const s = this.spectral();
    return { T: s.T, V: this.potentialEnergy(), E: s.T + this.potentialEnergy() };
  }

  // probability in a region (1D): Σ|ψ|²dx over x in [a, b]
  probIn(a, b) {
    let s = 0;
    for (let i = 0; i < this.n; i++) {
      const x = this.x(i);
      if (x >= a && x <= b) s += this.re[i] * this.re[i] + this.im[i] * this.im[i];
    }
    return s * this.dx;
  }

  // absorbing boundaries: cos²-ramped damping over the outer `frac` of the box
  setAbsorber(frac, strength) {
    this.mask = new Float64Array(this.size);
    const ramp = u => u <= 0 ? 0 : Math.sin(Math.min(1, u) * Math.PI / 2) ** 2;
    if (this.dim === 1) {
      const w = this.L * frac;
      for (let i = 0; i < this.n; i++) {
        const x = this.x(i);
        const e = Math.max(ramp((x - (this.L / 2 - w)) / w), ramp(((-this.L / 2 + w) - x) / w));
        this.mask[i] = Math.exp(-strength * e);
      }
    } else {
      const wx = this.Lx * frac, wy = this.Ly * frac;
      for (let j = 0; j < this.ny; j++) {
        for (let i = 0; i < this.nx; i++) {
          const x = this.x(i), y = this.y(j);
          const e = Math.max(
            ramp((x - (this.Lx / 2 - wx)) / wx), ramp(((-this.Lx / 2 + wx) - x) / wx),
            ramp((y - (this.Ly / 2 - wy)) / wy), ramp(((-this.Ly / 2 + wy) - y) / wy));
          this.mask[j * this.nx + i] = Math.exp(-strength * e);
        }
      }
    }
  }

  // two-particle (2D-as-x₁x₂) entanglement: purity of the reduced density
  // matrix ρ₁ = Tr₂|ψ⟩⟨ψ|;  effective Schmidt number K = 1/Tr ρ₁²
  purity() {
    const { nx, ny } = this;
    // B[a][b] = Σᵢ ψ(xᵢ,y_a)ψ*(xᵢ,y_b) — reduced matrix for particle 2
    let sum = 0;
    const Br = new Float64Array(ny * ny), Bi = new Float64Array(ny * ny);
    for (let a = 0; a < ny; a++) {
      for (let b = a; b < ny; b++) {
        let r = 0, m = 0;
        const oa = a * nx, ob = b * nx;
        for (let i = 0; i < nx; i++) {
          r += this.re[oa + i] * this.re[ob + i] + this.im[oa + i] * this.im[ob + i];
          m += this.re[oa + i] * this.im[ob + i] - this.im[oa + i] * this.re[ob + i];
        }
        Br[a * ny + b] = r; Bi[a * ny + b] = m;
        sum += (a === b ? 1 : 2) * (r * r + m * m);
      }
    }
    const norm = this.norm();
    // Tr ρ₂² = Σ|B|²·dx²·dy², with norm = Σ|ψ|²·dx·dy
    return sum * this.dx * this.dx * this.dy * this.dy / (norm * norm);
  }
}

// ---------------- eigenstates by imaginary-time relaxation ----------------

function solveEigen(state, count, opts = {}) {
  const dtau = opts.dtau || 2e-3;
  const maxIter = opts.maxIter || 6000;
  const tol = opts.tol || 1e-8;
  const found = [];
  const energies = [];
  for (let nEig = 0; nEig < count; nEig++) {
    const s = new QMState(state.dim === 1
      ? { dim: 1, n: state.n, L: state.L, m: state.m }
      : { dim: 2, nx: state.nx, ny: state.ny, Lx: state.Lx, Ly: state.Ly, m: state.m });
    s.V = state.V;
    s.oddProject = state.oddProject;
    // start from a randomized-but-deterministic seed with the right reach
    for (let i = 0; i < s.size; i++) {
      s.re[i] = Math.sin(1.7 * (i + 1) * (nEig + 1)) + 0.3 * Math.cos(0.37 * i);
      s.im[i] = 0;
    }
    if (s.oddProject) s.applyOdd();
    s.normalize();
    let prevE = Infinity;
    for (let it = 0; it < maxIter; it++) {
      s.imagStep(dtau);
      if (s.oddProject) s.applyOdd();
      for (const f of found) s.subtractProjection(f);
      s.normalize();
      if (it % 40 === 0) {
        const E = s.energy().E;
        if (Math.abs(E - prevE) < tol * Math.max(1, Math.abs(E))) break;
        prevE = E;
      }
    }
    found.push(s);
    energies.push(s.energy().E);
  }
  return { energies, states: found };
}

// analytic transmission through a rectangular barrier (height V0, width a)
function barrierT(E, V0, a, m = 1) {
  if (E <= 0) return 0;
  if (Math.abs(E - V0) < 1e-12) E += 1e-9;
  if (E < V0) {
    const kap = Math.sqrt(2 * m * (V0 - E));
    const sh = Math.sinh(kap * a);
    return 1 / (1 + V0 * V0 * sh * sh / (4 * E * (V0 - E)));
  }
  const k2 = Math.sqrt(2 * m * (E - V0));
  const sn = Math.sin(k2 * a);
  return 1 / (1 + V0 * V0 * sn * sn / (4 * E * (E - V0)));
}

// packet-averaged transmission: ∫|φ(k)|² T(k) dk for a Gaussian packet
function packetT(k0, sigmaX, V0, a, m = 1) {
  const sk = 1 / (2 * sigmaX);     // Δk of a Gaussian with Δx = σ
  let num = 0, den = 0;
  for (let i = -120; i <= 120; i++) {
    const k = k0 + 3 * sk * i / 120;
    const u = (k - k0) / sk;
    const w = Math.exp(-u * u / 2);
    num += w * barrierT(k * k / (2 * m), V0, a, m);
    den += w;
  }
  return num / den;
}

if (typeof module !== "undefined") {
  module.exports = { fft, kOf, QMState, solveEigen, barrierT, packetT };
}
