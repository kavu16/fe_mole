//! Interaction potentials and the force kernels that evaluate them.
//!
//! Kernels take plain slices rather than `&mut System`, so each one is
//! testable in isolation and parallelisable at M5 without touching storage.
//! M4 adds real-space Ewald and PME alongside [`lj`].

pub mod lj;
