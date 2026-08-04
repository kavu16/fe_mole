# fe_mole

A molecular dynamics engine written from scratch in Rust — Lennard-Jones and
electrostatics, periodic boundaries, validated against published results for
liquid argon.

**Status: M0 (skeleton and infrastructure).** There is no physics in the engine
yet. See [`MILESTONES.md`](MILESTONES.md) for the plan and the numeric
acceptance criteria that define "done" for each stage.

## What this is

A portfolio project, built to the standard that scientific software actually
requires: a physics change is not finished until a numerical check confirms it,
and no performance claim lands without a before/after measurement. Every
milestone has acceptance criteria stated as numbers rather than judgment calls,
and the measured values are recorded next to them.

## Build and test

```sh
cargo test                                  # unit and integration tests
cargo clippy --all-targets -- -D warnings   # lints, warnings are errors
cargo fmt --check                           # formatting
cargo bench                                 # criterion benchmarks
```

The toolchain is pinned in `rust-toolchain.toml`, so CI and a local checkout
build with the same compiler.

## Units

The engine core works in Å, amu, fs, kcal/mol, K and elementary charge
throughout. Conversion to any other system happens at the boundary. The
constants are in `src/units.rs` and each is checked against SI by a test that
derives it — see [`docs/theory/units.md`](docs/theory/units.md).

## Layout

```
src/units.rs       unit system as named constants, SI-checked
src/geometry.rs    Vec3, periodic box, minimum image convention
src/system.rs      particle storage (struct-of-arrays)
benches/           criterion benchmarks
docs/theory/       physics and numerics notes, in the code's notation
docs/decisions/    architecture decision records
validation/        plots and data backing each acceptance criterion
LOG.md             lab notebook: what was tried, what happened, the numbers
```

## Reference system

Validation targets Rahman's 1964 liquid argon simulation: `N = 864`,
`T = 94.4 K`, `ρ = 1.374 g/cm³`, `σ = 3.4 Å`, `ε/k_B = 120 K`.

A. Rahman, "Correlations in the Motion of Atoms in Liquid Argon",
*Phys. Rev.* **136**, A405 (1964).

## License

Not yet chosen.
