# CLAUDE.md

## What this project is

A molecular dynamics engine written from scratch in Rust. It is a learning
project: the goal is to understand the physics and numerics by implementing
them, not to compete with LAMMPS or GROMACS. That shapes every tradeoff:

- **Correctness beats speed, and validated correctness beats claimed correctness.**
  A physics change is not done until a numerical check confirms it. An engine
  you cannot verify teaches you nothing.
- **Clarity beats cleverness.** This code is meant to be re-read and understood
  months later. Prefer an obvious implementation with a comment citing the
  paper over a terse one.
- **Optimize only with measurements.** No performance change lands without a
  before/after benchmark number.

## How we work

Because this is a learning project, work is split by **learning value, not by
difficulty**. Implementing a Lennard-Jones kernel is the point; wiring up a
criterion benchmark is not. Do not optimise for finishing a milestone quickly.

**Default division of labour.** For each chunk of work, Claude writes:

- the failing acceptance test,
- function signatures with doc comments stating the physics, the units, and the
  source (paper or `docs/theory/` note),
- the skeleton of any new `docs/theory/` note,
- all scaffolding: CI, criterion boilerplate, plotting scripts, file I/O,
  manifest edits.

**I write the function bodies** — the physics and the numerics. Claude does not
write them, even when asked indirectly ("what would that look like?", "can you
show me roughly?"). Answer with the hint ladder below instead.

**Tests must assert physics, not an implementation.** A test that encodes
Claude's version of the algorithm teaches nothing and hides bugs the two share.
Prefer, in rough order of value:

- comparison against an analytical result (Madelung constant, virial theorem),
- a conservation law (energy, momentum, time-reversibility),
- a numerical derivative: any force must match `-(V(r+h) − V(r−h)) / 2h`,
- comparison against a slower, obviously-correct implementation,
- a symmetry or limiting case (`r → ∞`, `N = 2`, zero temperature).

**Work in checkpoints, not milestones.** Break a milestone into pieces that
each end at something runnable and checkable, and stop there for review. M1 is
at least three: minimum-image-based pair loop, LJ potential and force, velocity
Verlet. Never implement the integrator and the force calculation in one go —
that is also a debugging rule (see Pitfalls), not just a pedagogical one.

### Hint ladder

When I say I am stuck, I will name a level. Give **that level and stop** — do
not volunteer a higher one, and do not append the answer "for reference". If I
do not name a level, default to 1.

1. **Direction** — which function, concept, or invariant is implicated. No
   equations.
2. **Reference** — the specific equation, paper section, or `docs/theory/` note
   to read.
3. **Pseudocode** — structure and order of operations, in prose or maths. No
   Rust.
4. **Code** — the implementation.

The same ladder applies when reviewing my code and finding a *physics or
numerics* error: report the symptom and where it would show up, at level 1, and
let me find it. Mechanical errors — typos, borrow-checker complaints, clippy
lints, formatting — are not a learning opportunity; just say what is wrong and
fix it.

### Bug-injection drills

At the end of a milestone, when I ask for one, introduce a single realistic bug
on a scratch branch and report **only the symptom**: a drift number, a shifted
`g(r)` peak, a heating trend, a failed conservation check. Realistic means the
kinds listed under Pitfalls — a dropped unit conversion, a cutoff applied to
the energy but not the force, a sign error, a neighbour list that misses pairs,
a thermostat that biases the sampled distribution.

Do not reveal the bug, the file, or the diff until I have committed to a
diagnosis. Then show it and we compare reasoning. Record the drill in `LOG.md`
if the diagnosis was wrong — a wrong diagnosis is exactly the kind of dead end
that log is for.

### Style

Write idiomatic Rust. Where a Rust idiom and a Fortran/C++ MD idiom conflict,
prefer the Rust one unless it costs measurable performance.

**Comments: err on the side of fewer.** A comment earns its place by saying
something the code cannot — a citation, a unit, a non-obvious invariant, a
measured number, a rejected alternative. Delete anything that restates the
line below it. Long explanations belong in `docs/theory/` or an ADR, with the
code pointing at them.

**Attributes are not decoration.** Both of these were over-applied in M0 and
have since been trimmed:

- `#[inline]` — only on a function small enough that the call overhead rivals
  the body *and* on a path measured to be hot. LLVM already inlines within a
  crate; the attribute mainly enables inlining *across* crates, which
  `lto = "fat"` in the release profile largely does anyway. Applying it
  untested contradicts "optimize only with measurements".
- `#[must_use]` — only where discarding the result is a likely bug: pure
  transforms a caller might mistake for in-place mutation (`SimBox::wrap`,
  `minimum_image`), and constructors. Not on plain getters; std marks those
  because it is a stability-guaranteed public API, which this is not.

`missing_docs` is on, so every public item needs a doc comment. One line is a
complete answer for an obvious getter.

### Error handling

The engine panics rather than returning `Result`, and that is deliberate. The
deciding question is **where the value came from**, not what is being checked:

- **Panic** on a contract violation — a bug in the calling code.
  `SimBox::cubic(-1.0)` is not a runtime condition to handle; no `?` anywhere up
  the stack makes it recoverable.
- **`Result`** for anything originating outside the process: config files,
  coordinate files, command-line arguments.

Every value reaching a constructor today comes from code in this repo, so there
is nothing a `Result` could mean. Two reasons that stays true for the core:

1. **There is no degraded mode.** A zero mass has no sane recovery — you cannot
   continue a trajectory with one slightly-wrong particle and report a
   slightly-wrong diffusion coefficient. A `Result` would push a decision onto
   a caller whose only option is to abort, and add `?` noise to every call site.
2. **Silent wrongness is the enemy.** A zero mass gives infinite acceleration,
   then NaN positions, which propagate for 10⁵ steps and yield a plausible
   energy trace. Panicking at construction is the loudest, earliest failure
   available. This is why `SimBox::new` rejects a zero-width box instead of
   letting `minimum_image` quietly return NaN.

Three tools, three jobs:

- `assert!` — validate at construction and at API boundaries. Cold paths only.
- `debug_assert!` — invariants checked in hot paths, where the cost would show
  up in a benchmark (e.g. the struct-of-arrays length invariant, and Newton's
  third law inside the pair loop later).
- `Result` — the I/O boundary, when it arrives. Use *parse, don't validate*:
  the parser returns `Result` and checks the numbers, then hands clean values
  to constructors that can go on asserting.

Document every panic in a `# Panics` doc section (Rust API guideline C-FAILURE).

This is idiomatic Rust generally, not a scientific-computing exception —
`ndarray` panics on shape mismatch and `nalgebra` on dimension mismatch, while
`Result`-heavy crates are the ones eating untrusted input. Numeric kernels lean
panic; I/O leans `Result`.

**A case where the answer depends on the source:** the minimum image convention
is only valid for `r_c ≤ L/2`, so that check must exist before the LJ cutoff
does. It is a panic if the cutoff is a constant in the code, and a `Result` if
it comes from a config file — same check, different answer.

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
