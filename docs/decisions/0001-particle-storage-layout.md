# 0001 — Hand-rolled struct-of-arrays for particle storage

**Status:** Accepted
**Date:** 2026-08-03
**Milestone:** M0

## Context

Particle state is the hottest data in the engine. Every force evaluation
streams positions; every integrator step touches positions, velocities and
forces. The layout chosen at M0 constrains the SIMD work at M5, and changing it
later means rewriting every kernel.

Two requirements were binding:

1. **Force kernels need positions immutably and forces mutably at the same
   time.** Any layout that cannot express that disjoint borrow forces either
   `unsafe`, or a wasteful copy, or index-based access that defeats
   bounds-check elimination.
2. **The layout has to remain replaceable at M5.** Over-aligning an array for
   AVX-512 or reordering particles into cell order must not require touching
   call sites.

## Options

**A. `Vec<Particle>` (array-of-structs).** Simplest to read and write. A pair
loop touching only `r` still pulls `v`, `f`, `mass`, `charge` and `kind` into
cache with it — roughly 3× the useful traffic — and it cannot be vectorised
over particles without a gather. Rejected on cache behaviour; this is the case
`CLAUDE.md` explicitly rules out.

**B. `soa-rs` derive macro.** Derives SoA storage from a struct definition,
with `Vec`-like `push`/`insert`/`remove` that keep every array consistent for
free. The borrow concern turned out **not** to apply: `Soa::slices_mut()`
returns a struct holding every field as a disjoint `&mut [T]`, which expresses
requirement 1 cleanly. This was a genuine contender, not a straw man.

**C. `Vec<[f64; 3]>` per property.** What LAMMPS (`double **x`) and GROMACS
(`rvec[]`) actually do at the top level. All three components of one particle
land in a single cache line, which suits scalar pair loops well. Poor fit for
wide SIMD over particles without a transpose — GROMACS in fact abandons it in
its inner kernels, packing into `x[4]y[4]z[4]` cluster blocks instead.

**D. Hand-rolled: one `Vec<f64>` per component per quantity.** Twelve parallel
arrays behind private fields and slice accessors.

## Decision

**D.** `soa-rs` was rejected on three points specific to this project:

- **Alignment control at M5.** `soa-rs` packs all fields into one allocation at
  computed offsets; a field's start is 8-byte aligned, and where it lands
  depends on `N`. Owning the `Vec`s means an individual array can be swapped
  for an over-aligned allocation behind the same `&[f64]` accessor, with zero
  call-site churn.
- **Fields have different lifetimes.** `mass`, `charge` and `kind` are static
  after setup; `f` is zeroed and rebuilt every step; `r` is read by everything.
  A single derived struct binds them all to one access pattern.
- **The layout is the thing being learned.** Memory layout is one of the two or
  three decisions that actually determine how an MD engine performs, and
  delegating it to a derive macro skips the part worth understanding. The
  hand-rolled version is about sixty lines of obvious code, and writing it
  forces the questions this record answers.

## Consequences

**Easy.** Each array is an independent, contiguous `f64` run. Field-level
borrow splitting inside a method yields disjoint slices with no `unsafe`.
Kernels take plain `&[f64]` / `&mut [f64]` as free functions, which makes them
testable in isolation and trivially `rayon`-parallelisable at M5.

**Hard.** Every operation that changes the particle count has to touch twelve
arrays by hand — `push` is verbose, and `remove`/`swap_remove`/spatial sorting
will be too when they arrive. This is the real cost of the decision. It is
mitigated by keeping the fields private and asserting the equal-length
invariant (`debug_assert_consistent`) so a missed array fails a test rather
than corrupting a run.

Each disjoint-borrow pattern needs its own accessor written by hand:
`split_for_forces` exists now, and M1 will want an equivalent for the
integrator. This is deliberate — one method per genuine access pattern, rather
than one bundle that hands out everything.

**Revisit if:** the array count grows past roughly twenty and `push` becomes a
liability, or profiling at M5 shows the flat layout losing to a cluster-blocked
one (`x[4]y[4]z[4]`) in the inner kernel. The second is a real possibility —
GROMACS made exactly that move — and the accessor boundary is what would make
it tractable.
