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

## 2026-08-04 — Repo made private; branch protection lost

**Context:** M0 infrastructure. Decided to keep the work private until it is
ready to share.

**What I did:** Switched `kavu16/fe_mole` from public to private.

**Result:** CI still runs on every push and is green (private repos get 2000
free Actions minutes/month; a run takes ~1 min, so the quota is not a
constraint at this rate). Branch protection was silently dropped: on the free
plan, both branch protection rules and rulesets require a paid plan for
private repositories. `GET /repos/kavu16/fe_mole/branches/main/protection`
now returns 403.

**Interpretation:** The CI *signal* survived; only the *enforcement* was lost.
A red run no longer blocks anything — it just reports. That is a meaningful
weakening for a project whose whole premise is that validation gates progress,
so it should not be left implicit: `MILESTONES.md` claimed `main` was
protected and has been corrected.

The practical substitute while private is a local `pre-push` hook running the
same three checks. That is weaker than server-side enforcement — it is
bypassable with `--no-verify` and does not survive a fresh clone — but it
catches the realistic failure, which is me pushing without running clippy.

**Next:** Re-enable protection on the `check` status when the repo goes public
(the command is in this entry's commit). Unaffected for M1: the acceptance
criteria are numbers from test runs, not from CI configuration.

---

## 2026-08-04 — Back to public; protection restored

**Context:** M0 infrastructure. Reversed the previous entry's decision after
reframing the project description.

**What I did:** Removed the language describing this as a portfolio project
aimed at particular employers, restating it as a learning project. Made the
repo public again and reapplied branch protection on `main`.

**Result:** Protection is enforced again: required status check `check`,
`strict` (branch must be current before merge), no required reviews,
`enforce_admins` false, no force-push, no branch deletion. The 403 from the
previous entry is gone.

**Interpretation:** The two settings that could have locked me out of my own
repo are both deliberately off. Required *reviews* would be unworkable solo —
GitHub does not let you approve your own pull request. `enforce_admins` stays
false so a direct push to `main` is still possible when warranted; the
protection documents the standard rather than making the repo unusable.

The reframing was worth doing on its own terms. The engineering rules it was
originally used to justify — validate before claiming, measure before
optimising — do not depend on who is reading. Motivating them by an audience
was the weaker argument; an engine you cannot verify teaches you nothing.

**Next:** M1. No infrastructure work outstanding.

---

## 2026-08-05 — M1 checkpoint 1: four bugs, and what each one looked like

**Context:** M1 checkpoint 1 — fcc lattice, Maxwell–Boltzmann velocities,
and the kinetic energy / temperature / momentum observables. Recording the
four things I got wrong, because the *symptoms* are the reusable part.

**What I did / Result:**

1. **fcc corner site used the bare cell index, not index × cell edge.**
   Symptom: particle count, density, and every-site-inside-the-box all passed;
   nearest-neighbour distance came out 0.155 Å against an expected 4.087 Å, and
   coordination numbers were 125–195 instead of 12. The lesson is in the
   *pattern* of failure: a uniformly wrong lattice would have been wrong by a
   clean factor. Getting the count right while the spacing is garbage means
   some sites are placed correctly and others are not — 0.155 Å was just the
   nearest accidental approach between the correct group and the crammed one.

2. **`sigma` missing the amu·Å²/fs² → kcal/mol conversion.** Caught by code
   review, not by a test, and that is the interesting part. The rescale-to-
   target-temperature step divides out any *global* error in the drawn
   variance, so the final velocities were going to be exactly right regardless.
   The bug was latent: correct output, incoherent code, and it would have
   detonated the first time anyone made the rescale optional.

3. **`remove_center_of_mass_momentum` subtracted `m·vᵢ` per particle** instead
   of the global `v_com = P/M`. A units check catches it immediately —
   amu·Å/fs cannot be subtracted from Å/fs. The awkward `*vx -= m * *vx`
   double-deref that prompted me to ask about syntax was a symptom of the wrong
   expression, not a Rust problem; the correct version has no `*vx` on the
   right-hand side at all.

4. **Rescale factor used `T_target/T_current` instead of its square root.**
   Symptom: temperature 0.81% low, kinetic energy 2.2% high. I asked whether
   this was floating point. It was not, and the magnitude is how you tell —
   double-precision round-off is ~1e-16 relative, so anything above ~1e-12 on a
   sum of a few thousand terms is structural. The residual was small only
   because the draw already put `T_current` near the target, making
   `rescale ≈ 1`, where a wrong function of a number near 1 is also near 1.

**Interpretation:** Three of the four were caught by tests; the one that was
not (#2) was the one a downstream normalisation would have masked forever.
That is the generalisable lesson from this checkpoint: **a normalisation step
downstream of a calculation hides errors in that calculation.** The fix was to
split `draw_maxwell_boltzmann` and `rescale_to_temperature` into separate
public functions so the raw draw's absolute scale is directly assertable
before anything normalises it.

Unit errors in this system are never subtle — a missing conversion is a factor
of 2390, not a few percent. So "off by 2%" positively rules out a units bug,
which is a useful thing to know while triaging.

**Next:** Checkpoint 2 — Lennard-Jones potential and force. The same masking
risk shows up again at M4: PME's reciprocal-space scale is partly fixed by the
self-energy correction, so a wrong prefactor can look right in the total energy
while the forces are wrong. Test the piece before the thing that normalises it.

---
