//! `fe_mole` — a molecular dynamics engine written from scratch.
//!
//! # Units
//!
//! The engine core works in a single internal unit system throughout:
//! Å, amu, fs, kcal/mol, K, and elementary charge. Conversion to and from any
//! other system (SI, reduced Lennard-Jones units) happens at the boundary —
//! no core routine ever sees reduced units. The constants that define the
//! system live in [`units`]; the derivation is in `docs/theory/units.md`.
//!
//! # Particle storage
//!
//! Particle data is struct-of-arrays: [`system::System`] holds one `Vec<f64>`
//! per component per quantity, not a `Vec<Particle>`. This is deliberate, for
//! SIMD and cache behaviour — see `docs/decisions/0001-particle-storage-layout.md`.

pub mod geometry;
pub mod init;
pub mod observables;
pub mod system;
pub mod units;
