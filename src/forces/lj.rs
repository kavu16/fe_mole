//! The Lennard-Jones 12-6 pair potential, with a shifted cutoff.
//!
//! See `docs/theory/lennard-jones.md`. This module is the scalar pair
//! interaction only — one pair, one separation. The O(N²) loop over pairs is
//! checkpoint 3.

use crate::geometry::SimBox;
use crate::system::{Slices3, SlicesMut3, System};
use crate::units::BOLTZMANN;

/// A Lennard-Jones interaction with a cutoff, for one pair of species.
///
/// Energies are kcal/mol, lengths Å. `energy_shift` is precomputed at
/// construction so the hot path never re-evaluates `U` at the cutoff.
///
/// # Cutoff
///
/// The potential is **energy-shifted**: `U(r) − U(r_c)` inside the cutoff,
/// zero outside. An unshifted truncation leaves `U` discontinuous at `r_c`,
/// and every pair crossing the boundary injects energy — which shows up in an
/// NVE run as monotonic heating rather than as random drift. Working out what
/// the shift does to the force is the point of `docs/theory/lennard-jones.md`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LennardJones {
    sigma: f64,
    epsilon: f64,
    r_cut: f64,
    energy_shift: f64,
}

/// `(σ/r)⁶`. Every LJ exponent is even, so `r` itself is never needed.
fn sigma_over_r_pow6(sigma: f64, r_squared: f64) -> f64 {
    let s2 = sigma * sigma / r_squared;
    s2 * s2 * s2
}

/// Unshifted 12-6 potential, kcal/mol. `new` needs it at the cutoff, where
/// [`LennardJones::energy`] would be circular.
fn bare_lj_potential(sigma: f64, epsilon: f64, r_squared: f64) -> f64 {
    let s6 = sigma_over_r_pow6(sigma, r_squared);
    4.0 * epsilon * (s6 * s6 - s6)
}

impl LennardJones {
    /// Constructs an interaction with length parameter `sigma` [Å], well depth
    /// `epsilon` [kcal/mol], and cutoff `r_cut` [Å].
    ///
    /// # Panics
    ///
    /// Panics if any argument is not finite and strictly positive, or if
    /// `r_cut <= sigma` — a cutoff inside the repulsive core would discard the
    /// entire attractive well and is far more likely a units mistake than an
    /// intentional choice.
    ///
    /// The other constraint on `r_cut`, that it must not exceed half the
    /// shortest box length, cannot be checked here because this type does not
    /// know the box. It belongs with the pair loop at checkpoint 3.
    #[must_use]
    pub fn new(sigma: f64, epsilon: f64, r_cut: f64) -> Self {
        assert!(sigma.is_finite() && sigma > 0.0, "positive");
        assert!(epsilon.is_finite() && epsilon > 0.0, "positive");
        assert!(r_cut.is_finite() && r_cut > 0.0, "positive");
        assert!(r_cut > sigma, "r_cut > sigma");
        Self {
            sigma,
            epsilon,
            r_cut,
            energy_shift: bare_lj_potential(sigma, epsilon, r_cut * r_cut),
        }
    }

    /// Argon, as used throughout the milestone list: `σ = 3.4 Å`,
    /// `ε/k_B = 120 K`, cutoff `2.5σ` (Rahman 1964).
    #[must_use]
    pub fn argon() -> Self {
        Self::new(3.4, 120.0 * BOLTZMANN, 2.5 * 3.4)
    }

    /// Length parameter, Å.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// Well depth, kcal/mol.
    pub fn epsilon(&self) -> f64 {
        self.epsilon
    }

    /// Cutoff radius, Å.
    pub fn r_cut(&self) -> f64 {
        self.r_cut
    }

    /// The constant subtracted from the potential inside the cutoff, kcal/mol.
    pub fn energy_shift(&self) -> f64 {
        self.energy_shift
    }

    /// Pair potential energy at squared separation `r_squared` [Å²], kcal/mol.
    ///
    /// Takes `r²` rather than `r` because that is what the pair loop has:
    /// [`crate::geometry::Vec3::norm_squared`] avoids a square root, and the
    /// cutoff comparison works just as well squared.
    ///
    /// Returns exactly zero at and beyond the cutoff.
    #[must_use]
    pub fn energy(&self, r_squared: f64) -> f64 {
        if r_squared >= self.r_cut() * self.r_cut() {
            return 0.0;
        }
        bare_lj_potential(self.sigma(), self.epsilon(), r_squared) - self.energy_shift()
    }

    /// `−U'(r) / r` at squared separation `r_squared` [Å²], in kcal/mol/Å².
    ///
    /// This is the quantity a pair loop actually wants: multiplying it by the
    /// displacement vector `r_ij` [Å] gives the force vector [kcal/mol/Å]
    /// directly, with no square root and no normalisation.
    ///
    /// Positive means repulsive — the force on `j` points away from `i`.
    /// Returns exactly zero at and beyond the cutoff.
    #[must_use]
    pub fn force_over_r(&self, r_squared: f64) -> f64 {
        if r_squared >= self.r_cut() * self.r_cut() {
            return 0.0;
        }

        let s6 = sigma_over_r_pow6(self.sigma(), r_squared);
        (48.0 / r_squared) * self.epsilon() * (s6 * s6 - s6 / 2.0)
    }
}

/// Adds Lennard-Jones forces from every pair into `forces`, and returns the
/// total LJ potential energy [kcal/mol].
///
/// **Accumulates.** Forces are added to whatever is already there, so a caller
/// combining several interactions zeroes once (see [`System::zero_forces`])
/// and then calls each kernel in turn.
///
/// O(N²): every pair is visited. M2 replaces this with neighbour lists, and
/// this path becomes the reference the fast path is validated against — so it
/// stays, and it stays obvious rather than clever.
///
/// Every pair displacement goes through [`SimBox::minimum_image`].
///
/// # Panics
///
/// Panics if `positions` and `forces` disagree in length, or if the cutoff
/// exceeds half the shortest box length. The minimum image convention is only
/// valid for `r_c ≤ L/2`; beyond it a particle can interact with two images of
/// the same neighbour, or with itself, and the energy silently stops being the
/// one the model defines. The cutoff is a compile-time constant here, so this
/// is a contract violation rather than a runtime condition — see `CLAUDE.md`
/// on error handling.
#[expect(unused_variables, reason = "stub: body is checkpoint 3 work")]
pub fn accumulate_forces(
    lj: &LennardJones,
    sim_box: &SimBox,
    positions: Slices3<'_>,
    forces: SlicesMut3<'_>,
) -> f64 {
    todo!("M1 checkpoint 3")
}

/// Zeroes the force arrays, evaluates the Lennard-Jones interaction over all
/// pairs, and returns the total potential energy [kcal/mol].
///
/// The convenience wrapper the integrator will call each step.
pub fn compute_forces(system: &mut System, lj: &LennardJones) -> f64 {
    system.zero_forces();
    let sim_box = *system.sim_box();
    let (positions, forces) = system.split_for_forces();
    accumulate_forces(lj, &sim_box, positions, forces)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const SIGMA: f64 = 3.4;
    const R_CUT: f64 = 2.5 * SIGMA;

    fn argon() -> LennardJones {
        LennardJones::argon()
    }

    /// Potential as a function of `r` rather than `r²`, for readability in
    /// tests that sweep separations.
    fn v(lj: &LennardJones, r: f64) -> f64 {
        lj.energy(r * r)
    }

    /// Analytic force magnitude `−U'(r)`, reconstructed from `force_over_r`.
    fn f(lj: &LennardJones, r: f64) -> f64 {
        lj.force_over_r(r * r) * r
    }

    #[test]
    fn epsilon_matches_the_reference_well_depth() {
        let lj = argon();
        assert_relative_eq!(lj.sigma(), SIGMA);
        assert_relative_eq!(lj.r_cut(), R_CUT);
        // ε/k_B = 120 K, so ε ≈ 0.2385 kcal/mol.
        assert_relative_eq!(lj.epsilon(), 120.0 * BOLTZMANN, max_relative = 1e-12);
        assert_relative_eq!(lj.epsilon(), 0.238_464_5, max_relative = 1e-6);
    }

    #[test]
    fn force_matches_the_numerical_derivative_of_the_energy() {
        // THE test for this checkpoint, and the M1 acceptance criterion for
        // the pair interaction: any force must match -(U(r+h) - U(r-h)) / 2h.
        // It catches sign errors, dropped chain-rule factors, and an energy and
        // a force that disagree about the cutoff — independently of whether the
        // potential itself is the one we intended.
        let lj = argon();
        let h = 1e-5;

        // Stay clear of the cutoff: the shifted potential is continuous there
        // but the force is not, so a central difference straddling r_cut is
        // meaningless.
        let mut r = 0.8 * SIGMA;
        while r < R_CUT - 10.0 * h {
            let numerical = -(v(&lj, r + h) - v(&lj, r - h)) / (2.0 * h);
            assert_relative_eq!(f(&lj, r), numerical, max_relative = 1e-6, epsilon = 1e-9);
            r += 0.01;
        }
    }

    #[test]
    fn the_minimum_sits_at_the_expected_depth_and_separation() {
        // U is minimised at r = 2^(1/6) σ, where it equals -ε (before the
        // shift) and the force vanishes. Three independent facts from one
        // point, and none of them depend on the cutoff treatment except
        // through the shift.
        let lj = argon();
        let r_min = 2.0_f64.powf(1.0 / 6.0) * SIGMA;

        assert_relative_eq!(f(&lj, r_min), 0.0, epsilon = 1e-12);
        assert_relative_eq!(
            v(&lj, r_min),
            -lj.epsilon() - lj.energy_shift(),
            max_relative = 1e-12
        );
        // Sanity: it really is a minimum.
        assert!(v(&lj, r_min) < v(&lj, r_min - 0.1));
        assert!(v(&lj, r_min) < v(&lj, r_min + 0.1));
    }

    #[test]
    fn the_unshifted_potential_crosses_zero_at_sigma() {
        // U(σ) = 0 by construction, so after shifting it is exactly -U(r_c).
        let lj = argon();
        assert_relative_eq!(v(&lj, SIGMA), -lj.energy_shift(), max_relative = 1e-12);
    }

    #[test]
    fn the_shift_makes_the_energy_continuous_at_the_cutoff() {
        // The whole reason the shift exists: without it there is a step at
        // r_c, and every pair crossing it injects energy, which appears as
        // monotonic heating in an NVE run.
        let lj = argon();
        let approaching = v(&lj, R_CUT - 1e-9);
        assert_relative_eq!(approaching, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn everything_vanishes_beyond_the_cutoff() {
        let lj = argon();
        for r in [R_CUT, R_CUT + 1e-6, R_CUT * 2.0, 1e3] {
            assert_eq!(v(&lj, r), 0.0, "energy at r = {r}");
            assert_eq!(lj.force_over_r(r * r), 0.0, "force at r = {r}");
        }
    }

    #[test]
    fn repulsive_inside_sigma_and_attractive_beyond_the_minimum() {
        let lj = argon();
        let r_min = 2.0_f64.powf(1.0 / 6.0) * SIGMA;

        // Short range: steeply repulsive, force pushes the pair apart.
        for r in [0.8 * SIGMA, 0.9 * SIGMA, 0.99 * SIGMA] {
            assert!(f(&lj, r) > 0.0, "expected repulsion at r = {r}");
        }
        // Between the minimum and the cutoff: attractive.
        for r in [r_min + 0.1, 0.5 * (r_min + R_CUT), R_CUT - 0.1] {
            assert!(f(&lj, r) < 0.0, "expected attraction at r = {r}");
        }
    }

    #[test]
    fn repulsion_grows_as_r_to_the_minus_twelve() {
        // Halving the separation deep inside the core multiplies the energy by
        // 2^12 = 4096, since the r^-12 term dominates. Checks the exponent
        // without restating the formula.
        //
        // The exact ratio is 4096·(x⁶ − 1)/(x⁶ − 64) for x = σ/r, so "deep"
        // means x⁶ ≫ 64. At 0.35σ the correction is still 13%; at 0.12σ it is
        // 1.9e-4. Physically absurd separations, but this tests a functional
        // form, not a configuration.
        let lj = LennardJones::new(SIGMA, 120.0 * BOLTZMANN, 100.0);
        let r = 0.12 * SIGMA;
        let ratio = (v(&lj, r) + lj.energy_shift()) / (v(&lj, 2.0 * r) + lj.energy_shift());
        assert_relative_eq!(ratio, 4096.0, max_relative = 1e-3);
    }

    #[test]
    fn attraction_decays_as_r_to_the_minus_six() {
        // Far outside σ the r^-6 term dominates, so doubling r divides the
        // magnitude by 2^6 = 64.
        let lj = LennardJones::new(SIGMA, 120.0 * BOLTZMANN, 1e4);
        let r = 6.0 * SIGMA;
        let ratio = (v(&lj, r) + lj.energy_shift()) / (v(&lj, 2.0 * r) + lj.energy_shift());
        assert_relative_eq!(ratio, 64.0, max_relative = 1e-2);
    }

    #[test]
    fn energy_scales_linearly_with_epsilon() {
        // ε sets the depth of the well and nothing else; doubling it doubles
        // the potential everywhere.
        let a = LennardJones::new(SIGMA, 0.1, 20.0);
        let b = LennardJones::new(SIGMA, 0.2, 20.0);
        for r in [0.9 * SIGMA, SIGMA, 1.5 * SIGMA, 3.0 * SIGMA] {
            assert_relative_eq!(v(&b, r), 2.0 * v(&a, r), max_relative = 1e-12);
            assert_relative_eq!(f(&b, r), 2.0 * f(&a, r), max_relative = 1e-12);
        }
    }

    #[test]
    #[should_panic(expected = "positive")]
    fn negative_sigma_is_rejected() {
        let _ = LennardJones::new(-1.0, 0.2, 8.5);
    }

    #[test]
    #[should_panic(expected = "r_cut")]
    fn a_cutoff_inside_the_core_is_rejected() {
        let _ = LennardJones::new(3.4, 0.2, 2.0);
    }

    // ---- checkpoint 3: the O(N²) pair loop ----

    use crate::geometry::Vec3;
    use crate::init::fcc_lattice;
    use rand::{RngExt, SeedableRng};
    use rand_pcg::Pcg64Mcg;

    const ARGON_MASS: f64 = 39.948;
    const ARGON_DENSITY: f64 = 1.374;

    /// An fcc lattice with every site nudged off its symmetry point, so forces
    /// are non-trivial but no pair is anywhere near overlapping.
    fn perturbed_lattice(cells: usize, jitter: f64, seed: u64) -> System {
        let n = 4 * cells * cells * cells;
        let sim_box = SimBox::from_density(n, ARGON_MASS, ARGON_DENSITY);
        let base = fcc_lattice(cells, sim_box, ARGON_MASS);
        let mut rng = Pcg64Mcg::seed_from_u64(seed);
        let mut s = System::with_capacity(n, sim_box);
        let r = base.positions();
        for i in 0..base.len() {
            let d = Vec3::new(
                rng.random_range(-jitter..jitter),
                rng.random_range(-jitter..jitter),
                rng.random_range(-jitter..jitter),
            );
            s.push(
                Vec3::new(r.x[i], r.y[i], r.z[i]) + d,
                Vec3::ZERO,
                ARGON_MASS,
                0.0,
                0,
            );
        }
        s
    }

    fn two_particles(sim_box: SimBox, a: Vec3, b: Vec3) -> System {
        let mut s = System::with_capacity(2, sim_box);
        s.push(a, Vec3::ZERO, ARGON_MASS, 0.0, 0);
        s.push(b, Vec3::ZERO, ARGON_MASS, 0.0, 0);
        s
    }

    /// Total energy by an independent double loop, with no force bookkeeping.
    /// The slower, obviously-correct implementation the kernel is checked
    /// against.
    fn reference_energy(lj: &LennardJones, system: &System) -> f64 {
        let r = system.positions();
        let sim_box = *system.sim_box();
        let mut total = 0.0;
        for i in 0..system.len() {
            for j in (i + 1)..system.len() {
                let a = Vec3::new(r.x[i], r.y[i], r.z[i]);
                let b = Vec3::new(r.x[j], r.y[j], r.z[j]);
                total += lj.energy(sim_box.minimum_image(b - a).norm_squared());
            }
        }
        total
    }

    fn net_force(system: &System) -> Vec3 {
        let f = system.forces();
        Vec3::new(f.x.iter().sum(), f.y.iter().sum(), f.z.iter().sum())
    }

    #[test]
    fn net_force_vanishes() {
        // Newton's third law, summed: every pair contributes +f and -f, so the
        // total must cancel to round-off. This is the checkpoint criterion, and
        // at M1 it is what makes momentum conserved over 10^5 steps.
        let lj = argon();
        let mut s = perturbed_lattice(3, 0.3, 1);
        compute_forces(&mut s, &lj);

        let p = net_force(&s);
        assert!(p.x.abs() < 1e-11, "net fx = {}", p.x);
        assert!(p.y.abs() < 1e-11, "net fy = {}", p.y);
        assert!(p.z.abs() < 1e-11, "net fz = {}", p.z);
    }

    #[test]
    fn a_perfect_fcc_lattice_has_zero_force_on_every_particle() {
        // fcc is a Bravais lattice, so every site is a centre of inversion:
        // for each neighbour at +d there is one at -d, and central forces
        // cancel exactly. True for any isotropic pair potential and any
        // isotropic cutoff, so this is an analytic check with no reference
        // value needed.
        let lj = argon();
        let mut s = perturbed_lattice(3, 0.0, 2);
        compute_forces(&mut s, &lj);

        let f = s.forces();
        for i in 0..s.len() {
            let mag = Vec3::new(f.x[i], f.y[i], f.z[i]).norm();
            assert!(mag < 1e-12, "particle {i} feels {mag} kcal/mol/Å");
        }
    }

    #[test]
    fn a_single_pair_obeys_newtons_third_law() {
        let lj = argon();
        let sim_box = SimBox::cubic(40.0);
        let sep = 4.0;
        let mut s = two_particles(
            sim_box,
            Vec3::new(10.0, 10.0, 10.0),
            Vec3::new(10.0 + sep, 10.0, 10.0),
        );
        let energy = compute_forces(&mut s, &lj);

        let f = s.forces();
        // Equal and opposite, along the separation, and matching the scalar
        // pair function exactly.
        let expected = lj.force_over_r(sep * sep) * sep;
        assert_relative_eq!(f.x[0], -expected, max_relative = 1e-12);
        assert_relative_eq!(f.x[1], expected, max_relative = 1e-12);
        assert_relative_eq!(f.y[0], 0.0, epsilon = 1e-15);
        assert_relative_eq!(f.z[0], 0.0, epsilon = 1e-15);
        assert_relative_eq!(energy, lj.energy(sep * sep), max_relative = 1e-12);
    }

    #[test]
    fn energy_matches_an_independent_pair_sum() {
        let lj = argon();
        let mut s = perturbed_lattice(3, 0.3, 3);
        let from_kernel = compute_forces(&mut s, &lj);
        assert_relative_eq!(from_kernel, reference_energy(&lj, &s), max_relative = 1e-12);
    }

    #[test]
    fn forces_match_the_numerical_gradient_of_the_total_energy() {
        // The many-body version of checkpoint 2's check, and the strongest test
        // here: it validates the loop structure, the minimum image, the cutoff
        // and the third-law bookkeeping at once, against an energy computed by
        // a completely separate routine.
        let lj = argon();
        let mut s = perturbed_lattice(2, 0.3, 4);
        compute_forces(&mut s, &lj);
        let analytic: Vec<Vec3> = (0..s.len())
            .map(|i| {
                let f = s.forces();
                Vec3::new(f.x[i], f.y[i], f.z[i])
            })
            .collect();

        let h = 1e-6;
        for i in (0..s.len()).step_by(7) {
            for axis in 0..3 {
                let shift = |s: &mut System, d: f64| {
                    let p = s.positions_mut();
                    match axis {
                        0 => p.x[i] += d,
                        1 => p.y[i] += d,
                        _ => p.z[i] += d,
                    }
                };
                shift(&mut s, h);
                let plus = reference_energy(&lj, &s);
                shift(&mut s, -2.0 * h);
                let minus = reference_energy(&lj, &s);
                shift(&mut s, h);

                let numerical = -(plus - minus) / (2.0 * h);
                let a = match axis {
                    0 => analytic[i].x,
                    1 => analytic[i].y,
                    _ => analytic[i].z,
                };
                assert_relative_eq!(a, numerical, max_relative = 1e-5, epsilon = 1e-7);
            }
        }
    }

    #[test]
    fn the_minimum_image_is_applied() {
        // Two particles either side of a face are 2 Å apart through the
        // boundary and 18 Å apart the other way. Without minimum image they
        // would be outside the cutoff and feel nothing.
        let lj = argon();
        let sim_box = SimBox::cubic(20.0);
        let mut s = two_particles(sim_box, Vec3::new(1.0, 5.0, 5.0), Vec3::new(19.0, 5.0, 5.0));
        let energy = compute_forces(&mut s, &lj);

        assert_relative_eq!(energy, lj.energy(4.0), max_relative = 1e-12);
        assert!(
            s.forces().x[0].abs() > 0.0,
            "no interaction across the boundary"
        );
        assert_relative_eq!(s.forces().x[0], -s.forces().x[1], max_relative = 1e-12);
    }

    #[test]
    fn pairs_beyond_the_cutoff_contribute_nothing() {
        let lj = argon();
        let sim_box = SimBox::cubic(60.0);
        let mut s = two_particles(
            sim_box,
            Vec3::new(5.0, 5.0, 5.0),
            Vec3::new(5.0 + lj.r_cut() + 1.0, 5.0, 5.0),
        );
        let energy = compute_forces(&mut s, &lj);
        assert_eq!(energy, 0.0);
        assert_eq!(s.forces().x, [0.0, 0.0]);
    }

    #[test]
    fn translating_the_system_leaves_forces_and_energy_unchanged() {
        // Only displacements enter the physics, so a rigid shift -- including
        // one that pushes particles outside the primary cell -- must change
        // nothing.
        let lj = argon();
        let mut a = perturbed_lattice(3, 0.3, 5);
        let energy_a = compute_forces(&mut a, &lj);

        let mut b = perturbed_lattice(3, 0.3, 5);
        {
            let p = b.positions_mut();
            for k in 0..p.x.len() {
                p.x[k] += 137.0;
                p.y[k] -= 42.5;
                p.z[k] += 8.25;
            }
        }
        let energy_b = compute_forces(&mut b, &lj);

        assert_relative_eq!(energy_a, energy_b, max_relative = 1e-10);
        for i in 0..a.len() {
            assert_relative_eq!(a.forces().x[i], b.forces().x[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn accumulate_adds_rather_than_overwrites() {
        // The contract that lets several interactions share one force array.
        let lj = argon();
        let mut s = perturbed_lattice(2, 0.3, 6);
        let once = compute_forces(&mut s, &lj);
        let first: Vec<f64> = s.forces().x.to_vec();

        let sim_box = *s.sim_box();
        let (positions, forces) = s.split_for_forces();
        let again = accumulate_forces(&lj, &sim_box, positions, forces);

        assert_relative_eq!(once, again, max_relative = 1e-12);
        let doubled: Vec<f64> = s.forces().x.to_vec();
        for (&twice, &once) in doubled.iter().zip(&first) {
            assert_relative_eq!(twice, 2.0 * once, max_relative = 1e-12);
        }
    }

    #[test]
    #[should_panic(expected = "minimum image")]
    fn a_cutoff_wider_than_half_the_box_is_rejected() {
        // r_c = 8.5 Å needs L ≥ 17 Å.
        let lj = argon();
        let mut s = two_particles(
            SimBox::cubic(16.0),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(5.0, 1.0, 1.0),
        );
        let _ = compute_forces(&mut s, &lj);
    }
}
