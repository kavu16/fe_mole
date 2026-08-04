# CLAUDE.md

## What this project is

A molecular dynamics engine written from scratch in Rust. It is a portfolio
project targeting scientific software engineering roles (D. E. Shaw Research,
Schrödinger). That shapes every tradeoff:

- **Correctness beats speed, and validated correctness beats claimed correctness.**
  A physics change is not done until a numerical check confirms it.
- **Clarity beats cleverness.** This code will be read by hiring managers.
  Prefer an obvious implementation with a comment citing the paper over a
  terse one.
- **Optimize only with measurements.** No performance change lands without a
  before/after benchmark number.

## Where the context lives

Read these when relevant — do not assume their contents:

- `MILESTONES.md` — current milestone, and the numeric acceptance criteria that
  define "done" for each. Check here before proposing what to work on next.
- `LOG.md` — running lab notebook of what was tried and what happened.
  Append an entry after any significant debugging session or result.
- `docs/theory/` — my own distilled notes on the physics and numerics
  (Ewald/PME, integrators, thermostats). These use the same notation as the
  code. Prefer these over restating textbook derivations.
- `docs/decisions/` — architecture decision records. Add one when making a
  choice with real alternatives; keep it to half a page.

## Unit convention

**This is the single most common source of silent bugs in MD. Never guess.**

Internal units are:

| Quantity | Unit |
|---|---|
| Length | Å |
| Mass | amu |
| Time | fs |
| Energy | kcal/mol |
| Temperature | K |
| Charge | elementary charge (e) |

Derived constants in this system:

- Force: kcal/mol/Å
- `k_B = 0.0019872041 kcal/(mol·K)`
- Acceleration: `a [Å/fs²] = 4.184e-4 * F [kcal/mol/Å] / m [amu]`
  The `4.184e-4` conversion is required. Omitting it produces a simulation that
  looks plausible and is wrong.

Reduced (LJ) units may be used in tests that compare against reduced-unit
literature values, but conversion happens at the boundary — the engine core is
always in the table above.

## Conventions that are not obvious from the code

- **Particle data is struct-of-arrays**, not array-of-structs. This is
  deliberate, for SIMD and cache behavior. Do not "simplify" it to a
  `Vec<Particle>`.
- **Minimum image convention applies to every pair displacement.** There is a
  single helper for this; use it rather than open-coding the wrap.
- **Positions are never wrapped in place during integration** — wrapping breaks
  mean-squared-displacement and diffusion measurements. Wrap only at output.
- Neighbor lists use a skin distance; rebuild is triggered by accumulated
  displacement, not by a fixed step count.
- Prefer `f64` throughout. Do not introduce `f32` for speed without a
  documented accuracy check.

## Pitfalls specific to this codebase

- **Energy drift is the canary.** If total energy in an NVE run starts
  drifting after a change, stop and find the cause. Do not compensate with a
  thermostat.
- **Never change the integrator and the force calculation in the same commit.**
  It makes drift regressions impossible to bisect.
- Cutoffs need consistent treatment between energy and force, otherwise the
  system heats. If a shift or switching function is applied to the potential,
  the force must match.
- Random number streams must be seedable and reproducible; stochastic
  integrators are untestable otherwise.

## Workflow

- Run `cargo test`, `cargo clippy -- -D warnings`, and `cargo fmt --check`
  before proposing a commit.
- Run the NVE energy-conservation check before any commit that touches forces,
  the integrator, or boundary conditions.
- Benchmarks are `criterion`; report the delta, not just the new number.
- One branch per milestone. Tag commits where a validation target passes
  (e.g. `m2-rdf-validated`) so known-good states stay reachable.
- Commit messages: what changed and why, and the validation result if relevant.

## When suggesting physics or numerics

Cite the source (paper, or the relevant note in `docs/theory/`). If a method
has known variants, say which one is being implemented and why — e.g. smooth
PME (Essmann 1995) rather than the original PME (Darden 1993), BAOAB rather
than naive Langevin splitting. If unsure whether an approach is standard, say
so rather than presenting it as settled.
