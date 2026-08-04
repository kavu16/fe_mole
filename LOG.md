# Lab Notebook

Append-only. Newest entries at the bottom. Do not edit past entries — if
something turned out to be wrong, write a new entry saying so.

Write an entry when: a validation target passes or fails, a bug takes more
than an hour, a design decision gets made, a benchmark changes meaningfully,
or an assumption turns out to be wrong. Include the numbers. Include the dead
ends — the dead ends are the most useful part when writing this up later, and
they are what makes the record honest.

**Entry format:**

```
## YYYY-MM-DD — short title

**Context:** what I was working on, which milestone.
**What I did:** the change or experiment.
**Result:** the numbers. Actual measured values, not impressions.
**Interpretation:** what I think it means.
**Next:** what this implies for the next step.
```

---

## 2026-08-03 — M0 skeleton: storage layout, unit constants, benchmark harness

**Context:** M0. Nothing existed but a scratch `main.rs` that pushed one
particle into a `soa_rs::Soa` with `f32` fields. Goal was infrastructure that
makes every later milestone checkable: unit constants that can be trusted, a
CI gate, and a benchmark harness that produces a number.

**What I did:**

1. Evaluated `soa-rs` against a hand-rolled struct-of-arrays for particle
   storage, and dropped the dependency (ADR 0001).
2. Wrote `units.rs` with `BOLTZMANN` and `FORCE_TO_ACCEL`, plus tests that
   *derive* both from SI rather than restating the literal.
3. Wrote `geometry.rs` (`Vec3`, `SimBox`, the single minimum-image helper) and
   `system.rs` (the SoA store).
4. Benchmarked two forms of the minimum image convention.
5. Fixed two real defects in the existing tree: `criterion` was in
   `[dependencies]`, so it was being linked into the release binary; and
   `Particle` used `f32` against the `f64` rule in `CLAUDE.md`.

**Result:**

- 26 tests pass. `cargo clippy --all-targets -- -D warnings` and
  `cargo fmt --check` are clean.
- `FORCE_TO_ACCEL` reproduces from SI to better than `1e-12` relative. It is an
  exact identity, not a measurement: `N_A` cancels between "per mole" in the
  energy unit and the definition of the amu, and the thermochemical calorie is
  defined as exactly 4184 J, leaving `4184 × 10⁻⁷`.
- `BOLTZMANN` reproduces from SI to `8e-8` relative.
- Minimum image, 10⁴ displacements, aarch64 (M-series), `lto = "fat"`,
  `codegen-units = 1`:

  | form | time | per pair |
  |---|---|---|
  | `d - L·round(d/L)` (used) | 11.354 µs | **1.135 ns** |
  | branch form | 43.544 µs | 4.354 ns |

**Interpretation:** The `soa-rs` borrow problem I expected to be the deciding
factor turned out not to exist — `Soa::slices_mut()` hands back disjoint
per-field `&mut [T]`, which is exactly what a force kernel needs. The crate was
dropped on other grounds (alignment control for M5, fields with different
lifetimes, and wanting the layout to be a decision I can defend rather than one
I delegated). Worth recording that the obvious objection was wrong.

The benchmark result was the opposite of my expectation. I assumed the branch
form would be faster and that using `round` was paying for correctness — the
code comment said so before I measured. It is 3.8× *slower*: the conditionals
are unpredictable on real displacement data, while `round` compiles to a single
branchless instruction. There is no correctness-versus-speed tradeoff here at
all, and the comment claiming one has been corrected. This is a good reminder
that the "obviously cheaper" scalar form is often not, once the branch
predictor is involved.

1.135 ns/pair is now the floor on every force loop in the engine: M1's O(N²)
kernel over 864 particles does ~373k pair displacements, so ~0.42 ms/step of
irreducible minimum-image cost before any Lennard-Jones arithmetic.

**Next:** M1 — velocity Verlet and the LJ kernel, in separate commits so a
drift regression stays bisectable. The first thing that needs building is the
NVE energy-conservation harness, since it is the acceptance criterion for
everything after it.

---
