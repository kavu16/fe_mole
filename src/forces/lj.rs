//! The Lennard-Jones 12-6 pair potential, with a shifted cutoff.
//!
//! See `docs/theory/lennard-jones.md`. This module is the scalar pair
//! interaction only — one pair, one separation. The O(N²) loop over pairs is
//! checkpoint 3.

use crate::units::BOLTZMANN;

/// A Lennard-Jones interaction with a cutoff, for one pair of species.
///
/// Energies are kcal/mol, lengths Å. `energy_shift` is precomputed at
/// construction so the hot path never re-evaluates `V` at the cutoff.
///
/// # Cutoff
///
/// The potential is **energy-shifted**: `V(r) − V(r_c)` inside the cutoff,
/// zero outside. An unshifted truncation leaves `V` discontinuous at `r_c`,
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
    #[expect(unused_variables, reason = "stub: body is checkpoint 2 work")]
    #[must_use]
    pub fn new(sigma: f64, epsilon: f64, r_cut: f64) -> Self {
        todo!("M1 checkpoint 2")
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
    #[expect(unused_variables, reason = "stub: body is checkpoint 2 work")]
    #[must_use]
    pub fn energy(&self, r_squared: f64) -> f64 {
        todo!("M1 checkpoint 2")
    }

    /// `−V'(r) / r` at squared separation `r_squared` [Å²], in kcal/mol/Å².
    ///
    /// This is the quantity a pair loop actually wants: multiplying it by the
    /// displacement vector `r_ij` [Å] gives the force vector [kcal/mol/Å]
    /// directly, with no square root and no normalisation.
    ///
    /// Positive means repulsive — the force on `j` points away from `i`.
    /// Returns exactly zero at and beyond the cutoff.
    #[expect(unused_variables, reason = "stub: body is checkpoint 2 work")]
    #[must_use]
    pub fn force_over_r(&self, r_squared: f64) -> f64 {
        todo!("M1 checkpoint 2")
    }
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

    /// Analytic force magnitude `−V'(r)`, reconstructed from `force_over_r`.
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
        // the pair interaction: any force must match -(V(r+h) - V(r-h)) / 2h.
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
        // V is minimised at r = 2^(1/6) σ, where it equals -ε (before the
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
        // V(σ) = 0 by construction, so after shifting it is exactly -V(r_c).
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
        let lj = LennardJones::new(SIGMA, 120.0 * BOLTZMANN, 100.0);
        let r = 0.35 * SIGMA;
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
}
