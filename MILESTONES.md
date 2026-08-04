# Milestones

Each milestone is done when its **acceptance criteria** pass — not when the
code runs. Criteria are numeric on purpose: "done" should not be a judgment
call. Record the measured value next to each criterion when it passes, and tag
the commit.

Reference conditions used throughout (Rahman 1964, liquid argon):
`N = 864`, `T = 94.4 K`, `ρ = 1.374 g/cm³`, `σ = 3.4 Å`, `ε/k_B = 120 K`.

---

## M0 — Skeleton and infrastructure

- [x] Cargo workspace laid out; core types defined
      — single crate, `lib.rs` + thin `main.rs` driver (ADR 0002)
- [x] Particle storage is struct-of-arrays from the start
      — hand-rolled parallel `Vec<f64>`, `soa-rs` dropped (ADR 0001)
- [x] Unit system chosen, documented in `CLAUDE.md`, and encoded as named
      constants (not magic numbers) — `src/units.rs`, each constant derived
      from SI by a test rather than restated
- [x] CI running `cargo test`, `cargo clippy -- -D warnings`, `cargo fmt --check`
      — plus `cargo bench --no-run`. Branch protection is **not** enabled:
      the repo is private, and GitHub gates protection rules and rulesets on
      a paid plan for private repos. CI still runs on every push. Re-enable
      protection on the `check` status when the repo goes public.
- [x] `criterion` benchmark harness in place with one trivial benchmark
      — `benches/minimum_image.rs`
- [x] `LOG.md`, `docs/theory/`, `docs/decisions/` created

**Acceptance:** CI green on a trivial commit. Benchmark harness produces a
number. No physics yet.

*Measured:* CI green on the initial commit (run 30876592198). 26 tests pass;
clippy `-D warnings` and `fmt --check` clean. Benchmark: minimum image at
**1.135 ns/pair** over 10⁴ displacements (aarch64, `lto = "fat"`,
`codegen-units = 1`) — against 4.354 ns/pair for the branch form, i.e. the
correct version is 3.8× faster, not slower. `FORCE_TO_ACCEL` reproduces from
SI to `< 1e-12` relative; `BOLTZMANN` to `8e-8`.

---

## M1 — Minimal correct NVE integrator

Lennard-Jones only, velocity Verlet, periodic boundaries, O(N²) forces.

- [ ] Velocity Verlet implemented
- [ ] LJ potential and force, consistently cut off (shifted or switched)
- [ ] Minimum image convention, single shared helper
- [ ] Total energy logged every N steps

**Acceptance:**
- [ ] Relative total-energy drift `|ΔE| / |E₀| < 1e-4` over 10⁵ steps at
      `dt = 1 fs` (equivalently `dt* ≈ 0.005` reduced)
- [ ] Drift scales as `O(dt²)`: halving `dt` reduces drift by ≈4×
- [ ] Momentum conserved to machine precision (total `|p| < 1e-10` after 10⁵ steps)
- [ ] Energy-vs-time plot committed to `validation/`

*Measured: _______*

**Tag:** `m1-nve-conserving`

---

## M2 — Neighbor lists and first real validation

- [ ] Cell lists + Verlet neighbor lists with skin distance
- [ ] Displacement-triggered rebuild
- [ ] Radial distribution function `g(r)` computed
- [ ] Mean-squared displacement / diffusion coefficient computed

**Acceptance:**
- [ ] Forces from the neighbor-list path match the O(N²) path to `< 1e-10`
      relative, on a randomized configuration (this is the regression test that
      protects every later milestone)
- [ ] Runtime scales approximately linearly in N over N = 10³ → 10⁵
- [ ] `g(r)` for liquid argon reproduces Rahman 1964: first peak position
      within 2%, first peak height within 5%
- [ ] Self-diffusion coefficient within ~15% of the Rahman value
      (≈ 2.4e-5 cm²/s at the reference conditions)
- [ ] `g(r)` overlay against the reference committed to `validation/`

*Measured: _______*

**Tag:** `m2-rdf-validated`

---

## M3 — Thermostats

- [ ] Langevin thermostat, BAOAB splitting (Leimkuhler & Matthews)
- [ ] Nosé–Hoover (or Nosé–Hoover chain)
- [ ] Seedable, reproducible RNG streams
- [ ] Note in `docs/theory/` on why BAOAB over naive splitting

**Acceptance:**
- [ ] Mean temperature within 1% of target over a 10⁵-step production run
- [ ] Kinetic energy distribution matches the analytic Maxwell–Boltzmann
      result (chi-squared or K-S test passes at the 5% level)
- [ ] Temperature-control error for BAOAB scales better with `dt` than the
      naive scheme — demonstrate with a `dt` sweep
- [ ] `g(r)` from NVT agrees with `g(r)` from M2's NVE run within statistical
      error (ensembles must agree on structure)

*Measured: _______*

**Tag:** `m3-nvt-validated`

---

## M4 — Electrostatics / PME

Build up in order: bare Coulomb → Ewald → smooth PME (Essmann 1995).

- [ ] Direct-space Ewald sum with real/reciprocal splitting
- [ ] Smooth PME with B-spline interpolation and FFT (`rustfft`)
- [ ] Self-energy and excluded-pair corrections handled
- [ ] Note in `docs/theory/pme.md` matching code notation

**Acceptance:**
- [ ] Madelung constant for the NaCl lattice reproduced to `< 1e-5` relative
      (target: 1.747565) — this is an exact analytical check, no excuses
- [ ] PME energy matches converged direct Ewald to `< 1e-5` relative on a
      disordered charged system
- [ ] PME forces match numerical derivatives of the PME energy to `< 1e-6`
- [ ] Energy conservation preserved: NVE drift with PME still meets the M1 bound
- [ ] Cost scales as `O(N log N)`, demonstrated over N = 10³ → 10⁵

*Measured: _______*

**Tag:** `m4-pme-validated`

---

## M5 — Performance engineering

- [ ] Profiled (`samply` / `perf`); hot paths identified and documented
- [ ] `rayon` parallelism over force computation
- [ ] SIMD in the inner force loop
- [ ] Strong and weak scaling studies

**Acceptance:**
- [ ] All M1–M4 validation criteria still pass unchanged (correctness is not
      negotiable for speed)
- [ ] Strong scaling: ≥ 60% parallel efficiency at 8 threads
- [ ] Weak scaling curve produced and committed
- [ ] Benchmarked against LAMMPS or GROMACS on an identical system; the ratio
      is reported honestly, with an explanation of where the gap comes from
- [ ] Optimization history documented — what was tried, what the numbers were,
      including changes that did not work

*Measured: _______*

**Tag:** `m5-optimized`

---

## M6 — Writeup

- [ ] README with a clear description, build instructions, and a results summary
- [ ] Technical writeup: physics, numerical methods, every validation result,
      performance analysis
- [ ] All validation plots reproducible from committed scripts
- [ ] Honest limitations section — what this engine does not do

**Acceptance:** someone with a computational chemistry background can read the
writeup and tell that the engine is correct, without running it.

---

## Scope discipline

If time gets short, **M1, M2, and M4 are the non-negotiable core** — a
validated, PME-capable engine. M3 and M5 are what elevate it. Resist adding
features (bonded forces, barostats, new file formats) until the core is done
and validated; an engine that does one thing verifiably right is worth more
here than one that does five things unverified.
