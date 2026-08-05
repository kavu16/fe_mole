//! Geometry: the [`Vec3`] scalar value type and the periodic [`SimBox`].

use std::ops::{Add, Mul, Neg, Sub};

use crate::units::AVOGADRO;

/// A three-component vector in internal units (Å, or Å/fs, or kcal/mol/Å —
/// the type does not track which).
///
/// This is a **scalar value type**: single displacements, single positions,
/// values crossing an API boundary. Bulk particle data is stored
/// struct-of-arrays in [`crate::system::System`]. Do not build a `Vec<Vec3>`
/// of particle state — that reintroduces the array-of-structs layout this
/// engine deliberately avoids.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Self = Self {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// Constructs a vector from its components.
    #[inline]
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Squared Euclidean norm. Prefer this over [`Vec3::norm`] in pair loops:
    /// cutoff comparisons work just as well on `r²` and skip the square root.
    #[inline]
    #[must_use]
    pub fn norm_squared(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm(self) -> f64 {
        self.norm_squared().sqrt()
    }
}

impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Vec3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Add for Vec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Vec3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Vec3::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        rhs * self
    }
}

impl Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

/// An orthorhombic periodic simulation box, side lengths in Å.
///
/// Orthorhombic only for now. Triclinic cells need a 3×3 matrix and a
/// different minimum-image routine; that is deferred until something in the
/// milestone list actually requires it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimBox {
    lengths: Vec3,
}

impl SimBox {
    /// Constructs a box with the given side lengths in Å.
    ///
    /// # Panics
    ///
    /// Panics if any side length is not finite and strictly positive. A
    /// zero-width box would make [`SimBox::minimum_image`] produce `NaN`
    /// silently, which is far harder to diagnose later than a panic here.
    #[must_use]
    pub fn new(lx: f64, ly: f64, lz: f64) -> Self {
        assert!(
            lx.is_finite() && ly.is_finite() && lz.is_finite(),
            "box side lengths must be finite, got ({lx}, {ly}, {lz})"
        );
        assert!(
            lx > 0.0 && ly > 0.0 && lz > 0.0,
            "box side lengths must be positive, got ({lx}, {ly}, {lz})"
        );
        Self {
            lengths: Vec3::new(lx, ly, lz),
        }
    }

    /// Constructs a cubic box of side `l` Å.
    ///
    /// # Panics
    ///
    /// Panics if `l` is not finite and strictly positive.
    #[must_use]
    pub fn cubic(l: f64) -> Self {
        Self::new(l, l, l)
    }

    /// Constructs a cubic box sized to hold `n` particles of mass `mass` [amu]
    /// at mass density `density` [g/cm³].
    ///
    /// The reference conditions used throughout the milestone list
    /// (Rahman 1964, liquid argon: `n = 864`, `mass = 39.948`,
    /// `density = 1.374`) give `L ≈ 34.68 Å`.
    ///
    /// # Panics
    ///
    /// Panics if `n` is zero, or if `mass` or `density` is not finite and
    /// strictly positive.
    #[must_use]
    pub fn from_density(n: usize, mass: f64, density: f64) -> Self {
        assert!(n > 0, "cannot size a box for zero particles");
        assert!(
            mass.is_finite() && mass > 0.0,
            "mass must be finite and positive, got {mass}"
        );
        assert!(
            density.is_finite() && density > 0.0,
            "density must be finite and positive, got {density}"
        );

        // Total mass in grams: n particles of `mass` amu is n·mass g/mol.
        let grams = n as f64 * mass / AVOGADRO;
        // Volume in cm³, then in Å³ (1 cm³ = 1e24 Å³).
        let volume = grams / density * 1e24;
        Self::cubic(volume.cbrt())
    }

    /// The box side lengths, Å.
    pub fn lengths(&self) -> Vec3 {
        self.lengths
    }

    /// Box volume, Å³.
    #[must_use]
    pub fn volume(&self) -> f64 {
        self.lengths.x * self.lengths.y * self.lengths.z
    }

    /// Applies the minimum image convention to a displacement.
    ///
    /// **This is the single shared helper.** Every pair displacement in the
    /// engine goes through here; do not open-code the wrap anywhere else.
    ///
    /// Uses the rounding form `d -= L · round(d / L)`, which is correct for a
    /// displacement of any magnitude. The branch form
    /// (`if d > L/2 { d -= L }`) is only valid for `|d| < 3L/2` and fails
    /// silently outside it — positions in this engine are never wrapped in
    /// place, so displacements really can exceed that bound.
    ///
    /// The branch form is also *slower*, which was not the expected result:
    /// `benches/minimum_image.rs` measures 1.14 ns/pair for this version
    /// against 4.35 ns/pair for the branch version (aarch64, M0 baseline).
    /// `round` compiles to a branchless instruction, while the conditionals
    /// are unpredictable on real displacement data. There is no
    /// correctness-versus-speed tradeoff here to revisit later.
    ///
    /// See Allen & Tildesley, *Computer Simulation of Liquids*, §1.6.3.
    #[inline]
    #[must_use]
    pub fn minimum_image(&self, d: Vec3) -> Vec3 {
        Vec3 {
            x: d.x - self.lengths.x * (d.x / self.lengths.x).round(),
            y: d.y - self.lengths.y * (d.y / self.lengths.y).round(),
            z: d.z - self.lengths.z * (d.z / self.lengths.z).round(),
        }
    }

    /// Maps a position into the primary cell `[0, L)`.
    ///
    /// **Output only.** Never call this during integration: wrapping positions
    /// in place destroys mean-squared displacement and therefore the diffusion
    /// coefficient measured at M2.
    #[must_use]
    pub fn wrap(&self, r: Vec3) -> Vec3 {
        Vec3 {
            x: r.x.rem_euclid(self.lengths.x),
            y: r.y.rem_euclid(self.lengths.y),
            z: r.z.rem_euclid(self.lengths.z),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::{Rng, SeedableRng};
    use rand_pcg::Pcg64Mcg;

    fn seeded_rng() -> Pcg64Mcg {
        Pcg64Mcg::seed_from_u64(0x5EED_0001)
    }

    #[test]
    fn volume_of_a_cube() {
        assert_relative_eq!(SimBox::cubic(3.0).volume(), 27.0);
        assert_relative_eq!(SimBox::new(2.0, 3.0, 4.0).volume(), 24.0);
    }

    #[test]
    fn rahman_reference_box() {
        // n = 864 argon atoms at 1.374 g/cm³.
        let b = SimBox::from_density(864, 39.948, 1.374);
        assert_relative_eq!(b.lengths().x, 34.68, max_relative = 1e-3);

        // Round-trip the density back out of the box to check the conversion
        // in both directions.
        let grams = 864.0 * 39.948 / AVOGADRO;
        let density = grams / (b.volume() / 1e24);
        assert_relative_eq!(density, 1.374, max_relative = 1e-12);
    }

    #[test]
    fn minimum_image_leaves_short_displacements_untouched() {
        let b = SimBox::cubic(10.0);
        let d = Vec3::new(1.0, -2.5, 4.9);
        let m = b.minimum_image(d);
        // Bit-identical: no arithmetic should have been applied.
        assert_eq!(m, d);
    }

    #[test]
    fn minimum_image_wraps_long_displacements() {
        let b = SimBox::cubic(10.0);
        assert_relative_eq!(b.minimum_image(Vec3::new(6.0, 0.0, 0.0)).x, -4.0);
        assert_relative_eq!(b.minimum_image(Vec3::new(-6.0, 0.0, 0.0)).x, 4.0);
        // Many box widths away -- the branch form would fail here.
        assert_relative_eq!(b.minimum_image(Vec3::new(103.0, 0.0, 0.0)).x, 3.0);
        assert_relative_eq!(b.minimum_image(Vec3::new(-97.0, 0.0, 0.0)).x, 3.0);
    }

    #[test]
    fn minimum_image_result_is_always_in_half_box() {
        let b = SimBox::new(10.0, 13.5, 7.25);
        let l = b.lengths();
        let mut rng = seeded_rng();

        for _ in 0..100_000 {
            // Span many box widths in both directions.
            let d = Vec3::new(
                rng.random_range(-50.0 * l.x..50.0 * l.x),
                rng.random_range(-50.0 * l.y..50.0 * l.y),
                rng.random_range(-50.0 * l.z..50.0 * l.z),
            );
            let m = b.minimum_image(d);

            // Tolerance covers round-off in the multiply-subtract only.
            assert!(m.x.abs() <= l.x / 2.0 * (1.0 + 1e-12), "x: {d:?} -> {m:?}");
            assert!(m.y.abs() <= l.y / 2.0 * (1.0 + 1e-12), "y: {d:?} -> {m:?}");
            assert!(m.z.abs() <= l.z / 2.0 * (1.0 + 1e-12), "z: {d:?} -> {m:?}");
        }
    }

    #[test]
    fn minimum_image_differs_from_input_by_whole_box_widths() {
        // The defining property: the image displacement is congruent to the
        // original modulo L. This is what makes the convention physical.
        let b = SimBox::new(10.0, 13.5, 7.25);
        let l = b.lengths();
        let mut rng = seeded_rng();

        for _ in 0..100_000 {
            let d = Vec3::new(
                rng.random_range(-20.0 * l.x..20.0 * l.x),
                rng.random_range(-20.0 * l.y..20.0 * l.y),
                rng.random_range(-20.0 * l.z..20.0 * l.z),
            );
            let m = b.minimum_image(d);
            let n = (d.x - m.x) / l.x;
            assert_relative_eq!(n, n.round(), epsilon = 1e-9);
        }
    }

    #[test]
    fn minimum_image_is_idempotent() {
        let b = SimBox::new(10.0, 13.5, 7.25);
        let mut rng = seeded_rng();
        for _ in 0..10_000 {
            let d = Vec3::new(
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
            );
            let once = b.minimum_image(d);
            assert_eq!(b.minimum_image(once), once);
        }
    }

    #[test]
    fn minimum_image_is_antisymmetric() {
        // r_ij = -r_ji must hold, or Newton's third law breaks in the force
        // loop and momentum stops being conserved.
        let b = SimBox::new(10.0, 13.5, 7.25);
        let mut rng = seeded_rng();
        for _ in 0..10_000 {
            let d = Vec3::new(
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
            );
            let forward = b.minimum_image(d);
            let backward = b.minimum_image(Vec3::new(-d.x, -d.y, -d.z));
            // Exactly on the L/2 boundary both images are legitimate, so
            // compare magnitudes rather than requiring exact negation.
            assert_relative_eq!(forward.x.abs(), backward.x.abs(), epsilon = 1e-12);
            assert_relative_eq!(forward.y.abs(), backward.y.abs(), epsilon = 1e-12);
            assert_relative_eq!(forward.z.abs(), backward.z.abs(), epsilon = 1e-12);
        }
    }

    #[test]
    fn wrap_maps_into_primary_cell() {
        let b = SimBox::new(10.0, 13.5, 7.25);
        let l = b.lengths();
        let mut rng = seeded_rng();
        for _ in 0..10_000 {
            let r = Vec3::new(
                rng.random_range(-500.0..500.0),
                rng.random_range(-500.0..500.0),
                rng.random_range(-500.0..500.0),
            );
            let w = b.wrap(r);
            assert!((0.0..l.x).contains(&w.x), "{r:?} -> {w:?}");
            assert!((0.0..l.y).contains(&w.y), "{r:?} -> {w:?}");
            assert!((0.0..l.z).contains(&w.z), "{r:?} -> {w:?}");
        }
    }

    #[test]
    fn wrapping_does_not_change_pair_displacements() {
        // The reason wrapping is safe at output time: it does not affect any
        // minimum-image pair displacement.
        let b = SimBox::cubic(10.0);
        let mut rng = seeded_rng();
        for _ in 0..10_000 {
            let a = Vec3::new(
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
            );
            let c = Vec3::new(
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
            );
            let raw = b.minimum_image(Vec3::new(a.x - c.x, a.y - c.y, a.z - c.z));
            let (wa, wc) = (b.wrap(a), b.wrap(c));
            let wrapped = b.minimum_image(Vec3::new(wa.x - wc.x, wa.y - wc.y, wa.z - wc.z));
            assert_relative_eq!(raw.x.abs(), wrapped.x.abs(), epsilon = 1e-9);
            assert_relative_eq!(raw.y.abs(), wrapped.y.abs(), epsilon = 1e-9);
            assert_relative_eq!(raw.z.abs(), wrapped.z.abs(), epsilon = 1e-9);
        }
    }

    #[test]
    #[should_panic(expected = "must be positive")]
    fn zero_width_box_is_rejected() {
        let _ = SimBox::cubic(0.0);
    }

    #[test]
    #[should_panic(expected = "must be finite")]
    fn non_finite_box_is_rejected() {
        let _ = SimBox::cubic(f64::NAN);
    }

    #[test]
    fn norms() {
        let v = Vec3::new(3.0, 4.0, 12.0);
        assert_relative_eq!(v.norm_squared(), 169.0);
        assert_relative_eq!(v.norm(), 13.0);
        assert_relative_eq!(Vec3::ZERO.norm(), 0.0);
    }
}

/// M1 checkpoint 0: operator overloads for [`Vec3`].
///
/// These tests will not compile until `Add`, `Sub` and `Mul<f64>` are
/// implemented for `Vec3`. A compile failure is the intended starting state —
/// `cannot subtract `Vec3` from `Vec3`` means the red phase, not a mistake.
///
/// The properties asserted here are the ones the engine actually relies on,
/// not an arbitrary algebra checklist. In particular
/// `pair_displacement_matches_componentwise_form` is a differential test
/// against the existing hand-written form, so the operators can be adopted in
/// the M1 force loop without changing any number.
#[cfg(test)]
mod vec3_ops {
    use crate::geometry::{SimBox, Vec3};
    use approx::assert_relative_eq;
    use rand::{Rng, SeedableRng};
    use rand_pcg::Pcg64Mcg;

    /// Componentwise comparison, for properties that do not hold exactly in
    /// floating point.
    fn assert_close(a: Vec3, b: Vec3) {
        assert_relative_eq!(a.x, b.x, max_relative = 1e-12);
        assert_relative_eq!(a.y, b.y, max_relative = 1e-12);
        assert_relative_eq!(a.z, b.z, max_relative = 1e-12);
    }

    #[test]
    fn operators_are_componentwise() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(0.5, -1.0, 4.0);
        assert_eq!(a + b, Vec3::new(1.5, 1.0, 7.0));
        assert_eq!(a - b, Vec3::new(0.5, 3.0, -1.0));
        assert_eq!(a * 2.0, Vec3::new(2.0, 4.0, 6.0));
    }

    #[test]
    fn zero_is_the_additive_identity() {
        let a = Vec3::new(1.5, -2.5, 3.5);
        assert_eq!(a + Vec3::ZERO, a);
        assert_eq!(a - Vec3::ZERO, a);
        assert_eq!(a * 1.0, a);
        assert_eq!(a * 0.0, Vec3::ZERO);
        assert_eq!(a - a, Vec3::ZERO);
    }

    #[test]
    fn addition_commutes() {
        // Exact in IEEE 754: addition commutes even though it does not
        // associate.
        let a = Vec3::new(0.1, 0.2, 0.3);
        let b = Vec3::new(0.7, -0.4, 1e17);
        assert_eq!(a + b, b + a);
    }

    #[test]
    fn addition_associates_and_scaling_distributes() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-4.0, 5.0, -6.0);
        let c = Vec3::new(0.25, 0.5, 0.75);
        assert_close((a + b) + c, a + (b + c));
        assert_close((a + b) * 3.0, a * 3.0 + b * 3.0);
    }

    #[test]
    fn scaling_scales_the_norm() {
        let a = Vec3::new(3.0, 4.0, 12.0);
        for s in [0.0f64, 0.5, 1.0, -2.0, 1e6] {
            assert_relative_eq!((a * s).norm(), s.abs() * a.norm(), max_relative = 1e-12);
        }
    }

    #[test]
    fn displacement_is_antisymmetric() {
        // r_ij = -r_ji. Newton's third law depends on this holding exactly:
        // the M1 pair loop computes one displacement and applies +f to one
        // particle and -f to the other.
        let mut rng = Pcg64Mcg::seed_from_u64(0x5EED_0002);
        for _ in 0..10_000 {
            let a = Vec3::new(
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
            );
            let b = Vec3::new(
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
                rng.random_range(-100.0..100.0),
            );
            assert_eq!(b - a, (a - b) * -1.0);
        }
    }

    #[test]
    fn triangle_inequality() {
        let mut rng = Pcg64Mcg::seed_from_u64(0x5EED_0003);
        for _ in 0..10_000 {
            let a = Vec3::new(
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
            );
            let b = Vec3::new(
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
            );
            assert!((a + b).norm() <= a.norm() + b.norm() + 1e-12);
        }
    }

    #[test]
    fn pair_displacement_matches_componentwise_form() {
        // The differential test that makes adopting these operators safe: the
        // minimum-image pair displacement computed with `-` must equal the
        // hand-written componentwise form used everywhere today, bit for bit.
        let sim_box = SimBox::new(10.0, 13.5, 7.25);
        let mut rng = Pcg64Mcg::seed_from_u64(0x5EED_0004);
        for _ in 0..10_000 {
            let a = Vec3::new(
                rng.random_range(-50.0..50.0),
                rng.random_range(-50.0..50.0),
                rng.random_range(-50.0..50.0),
            );
            let b = Vec3::new(
                rng.random_range(-50.0..50.0),
                rng.random_range(-50.0..50.0),
                rng.random_range(-50.0..50.0),
            );
            let with_operator = sim_box.minimum_image(b - a);
            let componentwise = sim_box.minimum_image(Vec3::new(b.x - a.x, b.y - a.y, b.z - a.z));
            assert_eq!(with_operator, componentwise);
        }
    }
}
