//! Geometry: the [`Vec3`] scalar value type and the periodic [`SimBox`].

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
    /// Panics if `n` is zero, or if `mass` or `density` is not strictly
    /// positive.
    #[must_use]
    pub fn from_density(n: usize, mass: f64, density: f64) -> Self {
        assert!(n > 0, "cannot size a box for zero particles");
        assert!(mass > 0.0, "mass must be positive, got {mass}");
        assert!(density > 0.0, "density must be positive, got {density}");

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
