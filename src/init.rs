//! Initial configuration: lattice positions and thermal velocities.
//!
//! Nothing here is part of a trajectory — these run once, before step zero.

use rand::Rng;

use crate::geometry::SimBox;
use crate::system::System;

/// Builds a face-centred cubic lattice of `4 · cells_per_side³` particles
/// filling `sim_box`, every particle of mass `mass` [amu].
///
/// fcc is the equilibrium structure of solid argon and the conventional start
/// for a liquid-argon run: it reaches liquid density without overlapping
/// repulsive cores, so the first force evaluation is finite rather than
/// astronomical. The conventional cubic cell holds four atoms, at fractional
/// coordinates `(0,0,0)`, `(0,½,½)`, `(½,0,½)`, `(½,½,0)`.
///
/// The Rahman 1964 reference system is `cells_per_side = 6`, giving
/// `4 · 6³ = 864` particles — which is where that otherwise-arbitrary `N`
/// comes from.
///
/// The lattice must be commensurate with `sim_box`: the cell edge is
/// `L / cells_per_side` along each axis. If it is not, the periodic images
/// fail to line up and the structure has a seam at the boundary.
///
/// Velocities start at zero. Charge is `0.0` and kind is `0` for every
/// particle, since M1 is Lennard-Jones only.
///
/// # Panics
///
/// Panics if `cells_per_side` is zero, or `mass` is not finite and positive.
#[expect(unused_variables, reason = "stub: body is checkpoint 1 work")]
#[must_use]
pub fn fcc_lattice(cells_per_side: usize, sim_box: SimBox, mass: f64) -> System {
    todo!("M1 checkpoint 1")
}

/// Draws velocities from the Maxwell–Boltzmann distribution at `temperature`
/// [K], removes the net centre-of-mass momentum, and rescales so that the
/// instantaneous temperature equals `temperature`.
///
/// Each Cartesian velocity component is an independent Gaussian with zero mean
/// and variance `k_B T / m`. Mind the units — `k_B T / m` in
/// kcal/(mol·amu) is not Å²/fs²; the conversion is in [`crate::units`].
///
/// All three properties above must hold *simultaneously* on return. They are
/// not independent, so the order in which the three operations are applied
/// decides whether they do.
///
/// `rng` is borrowed rather than constructed here so the caller owns the
/// stream: reproducibility requires a seeded generator, and stochastic
/// integrators are untestable without it (see `CLAUDE.md`, Pitfalls).
///
/// # Panics
///
/// Panics if `temperature` is not finite and positive, or if the system has
/// fewer than two particles.
#[expect(unused_variables, reason = "stub: body is checkpoint 1 work")]
pub fn maxwell_boltzmann_velocities(system: &mut System, temperature: f64, rng: &mut impl Rng) {
    todo!("M1 checkpoint 1")
}

/// Subtracts the mass-weighted mean velocity from every particle, leaving zero
/// total momentum.
///
/// Total momentum is conserved by Newton's equations, so any drift present at
/// step zero persists for the entire run. Three things go wrong if it is left
/// in: it contributes to the measured kinetic energy without being thermal
/// motion, so the reported temperature is too high; it biases every observable
/// derived from that temperature; and it makes mean-squared displacement grow
/// as `t²` rather than `t`, destroying the diffusion coefficient at M2.
///
/// # Panics
///
/// Panics if the system is empty.
#[expect(unused_variables, reason = "stub: body is checkpoint 1 work")]
pub fn remove_center_of_mass_momentum(system: &mut System) {
    todo!("M1 checkpoint 1")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Vec3;
    use crate::observables::{temperature, total_momentum};
    use approx::assert_relative_eq;
    use rand::SeedableRng;
    use rand_pcg::Pcg64Mcg;

    // Rahman 1964 liquid argon.
    const N_CELLS: usize = 6;
    const N_PARTICLES: usize = 864;
    const ARGON_MASS: f64 = 39.948;
    const ARGON_DENSITY: f64 = 1.374;
    const ARGON_TEMPERATURE: f64 = 94.4;
    /// Argon's Lennard-Jones length parameter, Å.
    const ARGON_SIGMA: f64 = 3.4;

    fn reference_lattice() -> System {
        let sim_box = SimBox::from_density(N_PARTICLES, ARGON_MASS, ARGON_DENSITY);
        fcc_lattice(N_CELLS, sim_box, ARGON_MASS)
    }

    fn thermalised(seed: u64) -> System {
        let mut s = reference_lattice();
        maxwell_boltzmann_velocities(
            &mut s,
            ARGON_TEMPERATURE,
            &mut Pcg64Mcg::seed_from_u64(seed),
        );
        s
    }

    /// Number of neighbours of each particle within `cutoff`, under minimum
    /// image. O(N²) — fine for a test, and deliberately independent of the
    /// neighbour-list machinery that arrives at M2.
    fn neighbour_counts(system: &System, cutoff: f64) -> Vec<usize> {
        let r = system.positions();
        let sim_box = *system.sim_box();
        let n = system.len();
        let mut counts = vec![0usize; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let a = Vec3::new(r.x[i], r.y[i], r.z[i]);
                let b = Vec3::new(r.x[j], r.y[j], r.z[j]);
                if sim_box.minimum_image(b - a).norm() < cutoff {
                    counts[i] += 1;
                    counts[j] += 1;
                }
            }
        }
        counts
    }

    fn min_pair_distance(system: &System) -> f64 {
        let r = system.positions();
        let sim_box = *system.sim_box();
        let n = system.len();
        let mut min = f64::INFINITY;
        for i in 0..n {
            for j in (i + 1)..n {
                let a = Vec3::new(r.x[i], r.y[i], r.z[i]);
                let b = Vec3::new(r.x[j], r.y[j], r.z[j]);
                min = min.min(sim_box.minimum_image(b - a).norm());
            }
        }
        min
    }

    #[test]
    fn fcc_has_four_atoms_per_conventional_cell() {
        assert_eq!(reference_lattice().len(), N_PARTICLES);
        let b = SimBox::cubic(10.0);
        assert_eq!(fcc_lattice(1, b, 1.0).len(), 4);
        assert_eq!(fcc_lattice(3, b, 1.0).len(), 108);
    }

    #[test]
    fn fcc_reproduces_the_reference_density() {
        use crate::units::AVOGADRO;
        let s = reference_lattice();
        let grams = s.len() as f64 * ARGON_MASS / AVOGADRO;
        let density = grams / (s.sim_box().volume() / 1e24);
        assert_relative_eq!(density, ARGON_DENSITY, max_relative = 1e-12);
    }

    #[test]
    fn fcc_sites_lie_inside_the_box() {
        let s = reference_lattice();
        let l = s.sim_box().lengths();
        let r = s.positions();
        for i in 0..s.len() {
            assert!((0.0..l.x).contains(&r.x[i]), "particle {i} x = {}", r.x[i]);
            assert!((0.0..l.y).contains(&r.y[i]), "particle {i} y = {}", r.y[i]);
            assert!((0.0..l.z).contains(&r.z[i]), "particle {i} z = {}", r.z[i]);
        }
    }

    #[test]
    fn fcc_nearest_neighbour_distance_is_a_over_root_two() {
        // In fcc the closest pair is a corner atom and an adjacent face
        // centre, at a/√2 for conventional cell edge a.
        let s = reference_lattice();
        let a = s.sim_box().lengths().x / N_CELLS as f64;
        let expected = a / 2.0_f64.sqrt();
        assert_relative_eq!(min_pair_distance(&s), expected, max_relative = 1e-12);

        // ~4.09 Å at the reference density, comfortably outside argon's
        // repulsive core. This is why an fcc start does not blow up.
        assert!(
            expected > ARGON_SIGMA,
            "nearest neighbours at {expected} Å are inside σ = {ARGON_SIGMA} Å"
        );
    }

    #[test]
    fn fcc_coordination_number_is_twelve() {
        // The structural signature of fcc: 12 nearest neighbours, against 6
        // for simple cubic and 8 for bcc. Evaluated under minimum image, so a
        // lattice incommensurate with the box shows up here as particles near
        // the boundary having the wrong count.
        let s = reference_lattice();
        let a = s.sim_box().lengths().x / N_CELLS as f64;
        // Between the first shell (a/√2 ≈ 0.71a) and the second (a).
        let counts = neighbour_counts(&s, 0.85 * a);
        let wrong: Vec<_> = counts
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c != 12)
            .take(5)
            .collect();
        assert!(
            wrong.is_empty(),
            "particles without 12 neighbours: {wrong:?}"
        );
    }

    #[test]
    fn velocities_are_reproducible_from_a_seed() {
        let a = thermalised(7);
        let b = thermalised(7);
        assert_eq!(a.velocities().x, b.velocities().x);
        assert_eq!(a.velocities().y, b.velocities().y);
        assert_eq!(a.velocities().z, b.velocities().z);

        // ...and different seeds genuinely differ, so the above is not
        // passing because every velocity is zero.
        assert_ne!(a.velocities().x, thermalised(8).velocities().x);
    }

    #[test]
    fn momentum_is_zero_after_initialisation() {
        let s = thermalised(11);
        let p = total_momentum(&s);
        assert!(p.x.abs() < 1e-12, "px = {}", p.x);
        assert!(p.y.abs() < 1e-12, "py = {}", p.y);
        assert!(p.z.abs() < 1e-12, "pz = {}", p.z);
    }

    #[test]
    fn temperature_matches_the_target() {
        // This and `momentum_is_zero_after_initialisation` have to hold at the
        // same time; satisfying either one alone is easy.
        let s = thermalised(13);
        assert_relative_eq!(temperature(&s), ARGON_TEMPERATURE, max_relative = 1e-10);
    }

    #[test]
    fn removing_drift_leaves_relative_velocities_untouched() {
        // Subtracting a constant from every velocity cannot change any
        // difference between two of them.
        let mut s = thermalised(17);
        let before: Vec<f64> = (1..s.len())
            .map(|i| s.velocities().x[i] - s.velocities().x[0])
            .collect();
        remove_center_of_mass_momentum(&mut s);
        let after: Vec<f64> = (1..s.len())
            .map(|i| s.velocities().x[i] - s.velocities().x[0])
            .collect();
        for (b, a) in before.iter().zip(&after) {
            assert_relative_eq!(b, a, epsilon = 1e-15);
        }
    }

    #[test]
    fn removing_drift_is_idempotent() {
        let mut s = thermalised(19);
        remove_center_of_mass_momentum(&mut s);
        let p = total_momentum(&s);
        assert!(p.x.abs() < 1e-12 && p.y.abs() < 1e-12 && p.z.abs() < 1e-12);
    }

    #[test]
    fn equipartition_holds_across_species() {
        // ⟨v²⟩ must scale as 1/m at fixed temperature. With one species, an
        // error using m where 1/m belongs is hidden by the final rescaling;
        // two species with a 100× mass ratio expose it.
        let (light, heavy) = (1.0, 100.0);
        let n = 2000;
        let mut s = System::with_capacity(n, SimBox::cubic(20.0));
        for i in 0..n {
            let mass = if i % 2 == 0 { light } else { heavy };
            s.push(Vec3::ZERO, Vec3::ZERO, mass, 0.0, u16::from(i % 2 == 1));
        }
        maxwell_boltzmann_velocities(&mut s, 300.0, &mut Pcg64Mcg::seed_from_u64(23));

        let v = s.velocities();
        let mean_sq = |parity: usize| {
            let items: Vec<f64> = (0..n)
                .filter(|i| i % 2 == parity)
                .map(|i| v.x[i] * v.x[i] + v.y[i] * v.y[i] + v.z[i] * v.z[i])
                .collect();
            items.iter().sum::<f64>() / items.len() as f64
        };
        // Sampling error on each mean is ~2.6% at 1000 particles.
        assert_relative_eq!(mean_sq(0) / mean_sq(1), heavy / light, max_relative = 0.15);
    }

    #[test]
    fn velocity_components_are_isotropic_and_uncorrelated() {
        // No axis is special, and the components are independent draws. A
        // single draw reused across x, y and z would pass the temperature and
        // momentum tests but fails the correlation check here.
        let s = thermalised(29);
        let v = s.velocities();
        let n = s.len() as f64;
        let mean_sq = |c: &[f64]| c.iter().map(|x| x * x).sum::<f64>() / n;
        let (mxx, myy, mzz) = (mean_sq(v.x), mean_sq(v.y), mean_sq(v.z));

        // ~5% sampling error at N = 864; 20% is a loose smoke test. The
        // rigorous distributional check is M3's acceptance criterion.
        assert_relative_eq!(mxx / myy, 1.0, max_relative = 0.2);
        assert_relative_eq!(mxx / mzz, 1.0, max_relative = 0.2);

        let cross = |a: &[f64], b: &[f64]| {
            a.iter().zip(b).map(|(p, q)| p * q).sum::<f64>() / n / mxx.sqrt() / myy.sqrt()
        };
        assert!(
            cross(v.x, v.y).abs() < 0.15,
            "⟨vx·vy⟩ = {}",
            cross(v.x, v.y)
        );
        assert!(
            cross(v.x, v.z).abs() < 0.15,
            "⟨vx·vz⟩ = {}",
            cross(v.x, v.z)
        );
    }
}
