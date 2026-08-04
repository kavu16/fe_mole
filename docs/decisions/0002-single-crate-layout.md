# 0002 — Single crate rather than a Cargo workspace

**Status:** Accepted
**Date:** 2026-08-03
**Milestone:** M0

## Context

`MILESTONES.md` lists "Cargo workspace laid out" under M0. That wording is
loose, and the actual question is whether the engine should be split across
several crates now or kept as one.

## Options

**A. Multi-crate workspace** — `crates/fe-mole-core`, `crates/fe-mole-cli`,
perhaps `crates/fe-mole-io`. Enforces layering at the compiler level: the core
physically cannot depend on the CLI. Allows the core to be published
independently, and gives finer-grained incremental rebuilds.

**B. Single crate, `lib.rs` + a thin `main.rs`.** Modules give the same
layering by convention. Integration tests and benchmarks consume the library
through its public API exactly as an external user would.

## Decision

**B.** There is one consumer of the core — the driver binary — and no second
one on the milestone list through M6. A workspace would buy enforced layering
that module boundaries already provide at this size, at the cost of feature
unification quirks, an extra decision about where benchmarks live, and more
manifest surface for a reader to skim before reaching any physics.

The library/binary split does the part that actually matters: it forces the
engine to have a public API, and it means `benches/` and `tests/` exercise that
API rather than reaching into internals.

## Consequences

**Easy.** One manifest, one lockfile, one `cargo test`. Anything under `src/`
can be used from benchmarks and integration tests without further plumbing.

**Hard.** Nothing stops a module from reaching into another's internals via
`pub(crate)`; layering is a convention here, not a compiler guarantee. Worth
watching once the module count grows.

**Revisit if:** a second real consumer appears (a Python binding via `pyo3`, or
a separate analysis tool), or if compile times become a genuine irritation and
the dependency graph has a clean cut. Splitting later is mechanical — moving
modules between crates is far cheaper than the layout decision in
[0001](0001-particle-storage-layout.md).
