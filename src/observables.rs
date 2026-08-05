//! Instantaneous observables computed from particle state.

use crate::geometry::Vec3;
use crate::system::System;

/// Total kinetic energy, kcal/mol.
///
/// `KE = ½ Σᵢ mᵢ |vᵢ|²`. The sum is in amu·Å²/fs² and has to be converted —
/// see [`crate::units::MASS_VELOCITY_SQ_TO_ENERGY`]. Getting this wrong puts
/// kinetic and potential energy on different scales, and their sum, which M1
/// is judged on, becomes meaningless.
#[expect(unused_variables, reason = "stub: body is checkpoint 1 work")]
#[must_use]
pub fn kinetic_energy(system: &System) -> f64 {
    todo!("M1 checkpoint 1")
}

/// Instantaneous temperature, K.
///
/// From equipartition, `T = 2·KE / (N_dof · k_B)`.
///
/// `N_dof = 3N − 3`, not `3N`: total momentum is fixed at zero by
/// initialisation and conserved thereafter, so three degrees of freedom are
/// unavailable to thermal motion. Using `3N` reports a temperature low by
/// `(3N−3)/3N` — 0.35% at `N = 864`. Small enough to overlook, and large
/// enough to matter against M3's 1% acceptance criterion.
///
/// # Panics
///
/// Panics if the system has fewer than two particles, where `3N − 3` is not a
/// meaningful count of degrees of freedom.
#[expect(unused_variables, reason = "stub: body is checkpoint 1 work")]
#[must_use]
pub fn temperature(system: &System) -> f64 {
    todo!("M1 checkpoint 1")
}

/// Total momentum `Σᵢ mᵢ vᵢ`, amu·Å/fs.
///
/// Conserved exactly by Newton's third law, so this doubles as a check on the
/// force loop: M1 requires `|p| < 1e-10` after 10⁵ steps. A drift here after
/// the integrator lands means the pair loop is not applying equal and opposite
/// forces.
#[expect(unused_variables, reason = "stub: body is checkpoint 1 work")]
#[must_use]
pub fn total_momentum(system: &System) -> Vec3 {
    todo!("M1 checkpoint 1")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::SimBox;
    use crate::init::maxwell_boltzmann_velocities;
    use crate::units::BOLTZMANN;
    use approx::assert_relative_eq;
    use rand::SeedableRng;
    use rand_pcg::Pcg64Mcg;

    fn two_particles(m0: f64, v0: Vec3, m1: f64, v1: Vec3) -> System {
        let mut s = System::with_capacity(2, SimBox::cubic(10.0));
        s.push(Vec3::ZERO, v0, m0, 0.0, 0);
        s.push(Vec3::ZERO, v1, m1, 0.0, 0);
        s
    }

    #[test]
    fn momentum_of_a_known_system() {
        // Σ m v needs no unit conversion: amu × Å/fs is already the unit.
        let s = two_particles(
            2.0,
            Vec3::new(1.0, 0.0, -3.0),
            5.0,
            Vec3::new(0.0, 2.0, 1.0),
        );
        let p = total_momentum(&s);
        assert_relative_eq!(p.x, 2.0);
        assert_relative_eq!(p.y, 10.0);
        assert_relative_eq!(p.z, -1.0);
    }

    #[test]
    fn opposed_particles_carry_no_net_momentum() {
        let v = Vec3::new(1.5, -2.0, 0.25);
        let s = two_particles(4.0, v, 4.0, -v);
        let p = total_momentum(&s);
        assert_relative_eq!(p.x, 0.0);
        assert_relative_eq!(p.y, 0.0);
        assert_relative_eq!(p.z, 0.0);
    }

    #[test]
    fn a_system_at_rest_has_no_kinetic_energy() {
        let s = two_particles(1.0, Vec3::ZERO, 2.0, Vec3::ZERO);
        assert_relative_eq!(kinetic_energy(&s), 0.0);
    }

    #[test]
    fn kinetic_energy_is_quadratic_in_velocity() {
        // Pure scaling, so it holds whatever the unit conversion turns out to
        // be — doubling every velocity must quadruple the kinetic energy.
        let v0 = Vec3::new(1e-3, -2e-3, 0.5e-3);
        let v1 = Vec3::new(-0.5e-3, 1e-3, 2e-3);
        let single = kinetic_energy(&two_particles(39.948, v0, 12.0, v1));
        let doubled = kinetic_energy(&two_particles(39.948, v0 * 2.0, 12.0, v1 * 2.0));
        assert_relative_eq!(doubled, 4.0 * single, max_relative = 1e-12);
    }

    #[test]
    fn kinetic_energy_satisfies_equipartition() {
        // The check that pins the absolute scale, and therefore the unit
        // conversion: each of the 3N − 3 degrees of freedom carries ½k_B T.
        let n = 500;
        let mut s = System::with_capacity(n, SimBox::cubic(50.0));
        for _ in 0..n {
            s.push(Vec3::ZERO, Vec3::ZERO, 39.948, 0.0, 0);
        }
        let target = 120.0;
        maxwell_boltzmann_velocities(&mut s, target, &mut Pcg64Mcg::seed_from_u64(31));

        let dof = 3.0 * n as f64 - 3.0;
        let expected = 0.5 * dof * BOLTZMANN * target;
        assert_relative_eq!(kinetic_energy(&s), expected, max_relative = 1e-10);
    }

    #[test]
    fn temperature_is_consistent_with_kinetic_energy() {
        let n = 200;
        let mut s = System::with_capacity(n, SimBox::cubic(50.0));
        for _ in 0..n {
            s.push(Vec3::ZERO, Vec3::ZERO, 20.0, 0.0, 0);
        }
        maxwell_boltzmann_velocities(&mut s, 250.0, &mut Pcg64Mcg::seed_from_u64(37));

        let dof = 3.0 * n as f64 - 3.0;
        assert_relative_eq!(
            temperature(&s),
            2.0 * kinetic_energy(&s) / (dof * BOLTZMANN),
            max_relative = 1e-12
        );
    }

    #[test]
    fn temperature_is_quadratic_in_velocity() {
        let build = |scale: f64| {
            let mut s = System::with_capacity(3, SimBox::cubic(10.0));
            s.push(Vec3::ZERO, Vec3::new(1e-3, 0.0, 0.0) * scale, 10.0, 0.0, 0);
            s.push(Vec3::ZERO, Vec3::new(0.0, 2e-3, 0.0) * scale, 10.0, 0.0, 0);
            s.push(Vec3::ZERO, Vec3::new(0.0, 0.0, -1e-3) * scale, 10.0, 0.0, 0);
            s
        };
        assert_relative_eq!(
            temperature(&build(3.0)),
            9.0 * temperature(&build(1.0)),
            max_relative = 1e-12
        );
    }

    #[test]
    #[should_panic(expected = "degrees of freedom")]
    fn temperature_of_a_single_particle_is_rejected() {
        let mut s = System::with_capacity(1, SimBox::cubic(10.0));
        s.push(Vec3::ZERO, Vec3::new(1e-3, 0.0, 0.0), 1.0, 0.0, 0);
        let _ = temperature(&s);
    }
}
