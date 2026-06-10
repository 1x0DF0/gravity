# 🌌 Gravity Lab — a spacetime curvature instrument

An interactive N-body laboratory in **real physical units** (AU · day · M☉) that
makes *testable predictions* — and tests them, live, against analytic general
relativity and celestial mechanics. Zero dependencies, pure HTML/CSS/JS, runs
entirely offline.

![solar system](docs/preview.png)
![S2 around Sgr A*](docs/s2.png)

## Run it

Open `index.html` in any modern browser, or serve the folder:

```sh
python3 -m http.server 8000      # → http://localhost:8000
```

Run the validation suite (loads the actual sim code headlessly and checks it
against analytic physics):

```sh
node test/run.js
```

## The physics

Everything is integrated in astronomical units where G is the square of the
Gauss gravitational constant, so real ephemeris-style numbers come out:

| Model | Equation |
|---|---|
| Newtonian gravity | F = G·m₁·m₂ / rⁿ (n editable; Plummer softening ε optional) |
| 1PN relativity | pairwise EIH equations of motion (Will's standard form) |
| Black hole | Paczyński–Wiita Φ = −GM/(r−rₛ) — exact Schwarzschild ISCO at 6GM/c² |
| Light rays | null test rays, weak-field deflection (GR factor 2), captured inside 1.5rₛ |
| Time dilation | dτ/dt = √(1 + 2Φ/c² − v²/c²) |

**Integrators:** symplectic Euler (1st), velocity Verlet (2nd, symplectic),
Yoshida (4th, symplectic), classical RK4 (for velocity-dependent GR forces) —
selectable live, with adaptive substepping on close encounters.

**Measurements the instrument makes from its own data:**

- **Perihelion precession** — perihelion passages are detected with sub-step
  parabolic interpolation; the measured Δϖ per orbit is displayed next to the
  analytic 1PN prediction 6πGM/(c²a(1−e²)).
- **Empirical Kepler test** — orbital periods are measured by angle accumulation
  and plotted as T vs a (log–log) against the T = 2π√(a³/GM) line.
- **Conservation tracking** — relative drift of energy, angular momentum and
  linear momentum on a log scale, live. Switch integrators and watch the
  4th-order symplectic line drop by five decades.
- **Osculating elements** — a, e, T, speed, clock rate for any selected body.

## Scenarios

| Scenario | What it demonstrates |
|---|---|
| Solar system | Real masses, a, e for all 8 planets, started at perihelion |
| Mercury precession | 1PN advance vs the analytic prediction (c scaled ×0.02 to make ×2500 the real 0.104″/orbit visible; set c×1 for reality) |
| S2 around Sgr A* | The real 16-yr, e = 0.88 orbit at the Galactic Centre with real c — its ≈12′/orbit Schwarzschild precession is the effect GRAVITY measured in 2020 |
| Black hole ISCO | Paczyński–Wiita orbits: stable at 8GM/c², plunging inside 6GM/c², capture at 2rₛ |
| Sun–Jupiter Trojans | L1–L5 computed by root-finding; rotating-frame view shows tadpole libration |
| Binary + planet | Circumbinary dynamics |
| Three-body figure-8 | Chenciner–Montgomery choreography, solar masses at AU scale |
| Planetesimal accretion | Seeded RNG (reproducible), momentum-conserving merges |

## Views & tools

- **Top view** — spacetime grid pulled toward masses, heat-colored by potential.
- **Sheet view** — 3D embedding of Φ(x,y); drag to orbit the camera. Trails
  store the potential at the moment they were laid down and ride the sheet.
- **Rotating frame** — co-rotating frame of the two heaviest bodies; grid shows
  the effective potential Φ − ½ω²ρ², Lagrange points marked L1–L5. Launching a
  body in this view automatically adds the co-rotation velocity, so you can
  park a probe at L4 by just clicking there.
- **💡 Light rays** — emit a fan of null rays across the view; paths persist as
  fading ghosts. Around the black hole you get visible lensing and capture.
- **🚀 Launch mode** — drag = position + velocity (shown in km/s).
- **Body editor** — type exact m, x, y, vx, vy for any selected body.
- **Data** — CSV export of recorded trajectories, JSON save/load of exact state.

## Validation

`test/run.js` checks the shipped code against analytic results — among them:

- Mercury's 1PN perihelion precession matches 6πGM/(c²a(1−e²)) to < 3%
  (and vanishes with relativity off)
- photon deflection matches 4GM/(c²b) to < 5%
- PW circular orbits: stable at 8GM/c², unstable at 5.2GM/c² (ISCO at 6)
- L1 at one Hill radius, L4 equilateral, probe parked at L4 librates
- figure-8 choreography returns to its start after one analytic period
- Yoshida-4 energy error ≪ Verlet at equal Δt; ΔL/L ~ 10⁻¹⁵
- Earth's clock runs slow by ≈ 1.5×10⁻⁸, measured periods obey Kepler to < 1%

## Things to try

- Load **Mercury precession**, watch Δϖ measured converge to Δϖ predicted in
  the Selected-body panel, then set the c multiplier to 1 and see the real
  0.10″/orbit prediction appear.
- Load **Trojans**, enable 🚀 Launch, click exactly on L4 — the probe stays,
  librating. Click slightly off — tadpole orbit.
- Load **ISCO**, press 💡 — photon paths bend around the hole; the innermost
  rays spiral in and vanish.
- Set force exponent **n = 2.1** in the solar system: orbits precess (Bertrand's
  theorem) and the Kepler plot walks off the 3/2 line.
