//! Thin driver for the `fe_mole` engine.
//!
//! At M0 this only builds the reference system and reports its geometry —
//! there is no physics yet. It exists so the binary target compiles and so
//! there is somewhere obvious for the M1 integrator loop to land.

use fe_mole::geometry::{SimBox, Vec3};
use fe_mole::system::System;
use fe_mole::units;

/// Rahman 1964 liquid argon, the reference system used throughout the
/// milestone list.
const N_PARTICLES: usize = 864;
/// Argon mass, amu.
const ARGON_MASS: f64 = 39.948;
/// Reference density, g/cm³.
const ARGON_DENSITY: f64 = 1.374;
/// Reference temperature, K.
const ARGON_TEMPERATURE: f64 = 94.4;

fn main() {
    let sim_box = SimBox::from_density(N_PARTICLES, ARGON_MASS, ARGON_DENSITY);
    let mut system = System::with_capacity(N_PARTICLES, sim_box);

    // Placeholder configuration: all particles at the origin. M1 replaces this
    // with an fcc lattice and a Maxwell-Boltzmann velocity draw.
    for _ in 0..N_PARTICLES {
        system.push(Vec3::ZERO, Vec3::ZERO, ARGON_MASS, 0.0, 0);
    }

    let l = sim_box.lengths();
    println!("fe_mole — M0 skeleton (no physics yet)");
    println!("  particles      {}", system.len());
    println!("  box            {:.4} × {:.4} × {:.4} Å", l.x, l.y, l.z);
    println!("  volume         {:.2} Å³", sim_box.volume());
    println!("  density        {ARGON_DENSITY} g/cm³");
    println!(
        "  k_B·T          {:.6} kcal/mol  (T = {ARGON_TEMPERATURE} K)",
        units::BOLTZMANN * ARGON_TEMPERATURE
    );
}
