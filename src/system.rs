//! Particle storage, struct-of-arrays.

use crate::geometry::{SimBox, Vec3};

/// Immutable per-component slices of a vector quantity over all particles.
///
/// Each slice has length [`System::len`]; element `i` of each belongs to
/// particle `i`.
#[derive(Debug, Clone, Copy)]
pub struct Slices3<'a> {
    /// x components.
    pub x: &'a [f64],
    /// y components.
    pub y: &'a [f64],
    /// z components.
    pub z: &'a [f64],
}

/// Mutable per-component slices of a vector quantity over all particles.
///
/// The three slices are disjoint, which is what lets a force kernel take
/// positions immutably and forces mutably at the same time.
#[derive(Debug)]
pub struct SlicesMut3<'a> {
    /// x components.
    pub x: &'a mut [f64],
    /// y components.
    pub y: &'a mut [f64],
    /// z components.
    pub z: &'a mut [f64],
}

/// The particle system: positions, velocities, forces and per-particle
/// properties, plus the periodic box they live in.
///
/// Storage is struct-of-arrays — one `Vec<f64>` per component per quantity,
/// not a `Vec<Particle>`. Access goes through `&[f64]` accessors, so an
/// individual array can be swapped for an over-aligned allocation at M5
/// without touching call sites. See ADR 0001.
///
/// All twelve arrays share a length, [`System::len`]; fields are private so
/// that invariant cannot be broken from outside.
///
/// **Positions are unwrapped.** A particle that has crossed the box fifty
/// times has a coordinate fifty box widths out. Wrapping in place would
/// destroy mean-squared displacement and the diffusion coefficient measured at
/// M2 — use [`SimBox::wrap`] at output only.
#[derive(Debug, Clone)]
pub struct System {
    // Positions, Å. Unwrapped -- see the type-level note above.
    rx: Vec<f64>,
    ry: Vec<f64>,
    rz: Vec<f64>,
    // Velocities, Å/fs.
    vx: Vec<f64>,
    vy: Vec<f64>,
    vz: Vec<f64>,
    // Forces, kcal/mol/Å. Zeroed and rebuilt every step.
    fx: Vec<f64>,
    fy: Vec<f64>,
    fz: Vec<f64>,
    // Static after setup.
    mass: Vec<f64>,
    charge: Vec<f64>,
    kind: Vec<u16>,

    sim_box: SimBox,
}

impl System {
    /// Creates an empty system in `sim_box`, preallocated for `n` particles.
    #[must_use]
    pub fn with_capacity(n: usize, sim_box: SimBox) -> Self {
        Self {
            rx: Vec::with_capacity(n),
            ry: Vec::with_capacity(n),
            rz: Vec::with_capacity(n),
            vx: Vec::with_capacity(n),
            vy: Vec::with_capacity(n),
            vz: Vec::with_capacity(n),
            fx: Vec::with_capacity(n),
            fy: Vec::with_capacity(n),
            fz: Vec::with_capacity(n),
            mass: Vec::with_capacity(n),
            charge: Vec::with_capacity(n),
            kind: Vec::with_capacity(n),
            sim_box,
        }
    }

    /// Appends one particle. Force components start at zero.
    ///
    /// `r` is in Å, `v` in Å/fs, `mass` in amu, `charge` in elementary
    /// charges. `kind` indexes a species table (introduced at M1).
    ///
    /// # Panics
    ///
    /// Panics if `mass` is not strictly positive — a zero or negative mass
    /// gives an infinite or sign-flipped acceleration, which shows up much
    /// later as an unexplained energy blow-up.
    pub fn push(&mut self, r: Vec3, v: Vec3, mass: f64, charge: f64, kind: u16) {
        assert!(
            mass > 0.0 && mass.is_finite(),
            "particle mass must be finite and positive, got {mass}"
        );
        self.rx.push(r.x);
        self.ry.push(r.y);
        self.rz.push(r.z);
        self.vx.push(v.x);
        self.vy.push(v.y);
        self.vz.push(v.z);
        self.fx.push(0.0);
        self.fy.push(0.0);
        self.fz.push(0.0);
        self.mass.push(mass);
        self.charge.push(charge);
        self.kind.push(kind);
        self.debug_assert_consistent();
    }

    /// The number of particles.
    pub fn len(&self) -> usize {
        self.rx.len()
    }

    /// Whether the system holds no particles.
    pub fn is_empty(&self) -> bool {
        self.rx.is_empty()
    }

    /// The periodic box.
    pub fn sim_box(&self) -> &SimBox {
        &self.sim_box
    }

    /// Unwrapped positions, Å.
    pub fn positions(&self) -> Slices3<'_> {
        Slices3 {
            x: &self.rx,
            y: &self.ry,
            z: &self.rz,
        }
    }

    /// Velocities, Å/fs.
    pub fn velocities(&self) -> Slices3<'_> {
        Slices3 {
            x: &self.vx,
            y: &self.vy,
            z: &self.vz,
        }
    }

    /// Forces, kcal/mol/Å.
    pub fn forces(&self) -> Slices3<'_> {
        Slices3 {
            x: &self.fx,
            y: &self.fy,
            z: &self.fz,
        }
    }

    /// Particle masses, amu.
    pub fn masses(&self) -> &[f64] {
        &self.mass
    }

    /// Particle charges, elementary charges.
    pub fn charges(&self) -> &[f64] {
        &self.charge
    }

    /// Per-particle species indices.
    pub fn kinds(&self) -> &[u16] {
        &self.kind
    }

    /// Clears the force arrays. Call once at the top of every force
    /// evaluation, before any kernel accumulates into them.
    pub fn zero_forces(&mut self) {
        self.fx.fill(0.0);
        self.fy.fill(0.0);
        self.fz.fill(0.0);
    }

    /// Borrows positions immutably and forces mutably at once, so a force
    /// kernel can read one while accumulating into the other.
    ///
    /// Two separate calls would conflict; splitting inside one method lets
    /// field-level borrow splitting see the six arrays are disjoint. Later
    /// milestones add further split accessors as they need them.
    pub fn split_for_forces(&mut self) -> (Slices3<'_>, SlicesMut3<'_>) {
        (
            Slices3 {
                x: &self.rx,
                y: &self.ry,
                z: &self.rz,
            },
            SlicesMut3 {
                x: &mut self.fx,
                y: &mut self.fy,
                z: &mut self.fz,
            },
        )
    }

    /// Checks the equal-length invariant. Debug-only: it runs on every `push`,
    /// and private fields already guarantee it structurally.
    fn debug_assert_consistent(&self) {
        debug_assert!(
            [
                self.ry.len(),
                self.rz.len(),
                self.vx.len(),
                self.vy.len(),
                self.vz.len(),
                self.fx.len(),
                self.fy.len(),
                self.fz.len(),
                self.mass.len(),
                self.charge.len(),
                self.kind.len(),
            ]
            .iter()
            .all(|&n| n == self.rx.len()),
            "struct-of-arrays length invariant violated"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn two_particle_system() -> System {
        let mut s = System::with_capacity(2, SimBox::cubic(10.0));
        s.push(
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(0.1, 0.2, 0.3),
            39.948,
            0.0,
            0,
        );
        s.push(
            Vec3::new(4.0, 5.0, 6.0),
            Vec3::new(0.4, 0.5, 0.6),
            1.008,
            -1.0,
            1,
        );
        s
    }

    #[test]
    fn push_keeps_every_array_the_same_length() {
        let s = two_particle_system();
        assert_eq!(s.len(), 2);
        assert!(!s.is_empty());
        for slice in [s.positions(), s.velocities(), s.forces()] {
            assert_eq!(slice.x.len(), 2);
            assert_eq!(slice.y.len(), 2);
            assert_eq!(slice.z.len(), 2);
        }
        assert_eq!(s.masses().len(), 2);
        assert_eq!(s.charges().len(), 2);
        assert_eq!(s.kinds().len(), 2);
    }

    #[test]
    fn empty_system_is_empty() {
        let s = System::with_capacity(16, SimBox::cubic(10.0));
        assert_eq!(s.len(), 0);
        assert!(s.is_empty());
    }

    #[test]
    fn push_stores_values_componentwise() {
        let s = two_particle_system();
        assert_eq!(s.positions().x, [1.0, 4.0]);
        assert_eq!(s.positions().y, [2.0, 5.0]);
        assert_eq!(s.positions().z, [3.0, 6.0]);
        assert_eq!(s.velocities().x, [0.1, 0.4]);
        assert_eq!(s.masses(), [39.948, 1.008]);
        assert_eq!(s.charges(), [0.0, -1.0]);
        assert_eq!(s.kinds(), [0, 1]);
    }

    #[test]
    fn forces_start_at_zero() {
        let s = two_particle_system();
        assert_eq!(s.forces().x, [0.0, 0.0]);
        assert_eq!(s.forces().y, [0.0, 0.0]);
        assert_eq!(s.forces().z, [0.0, 0.0]);
    }

    #[test]
    fn zero_forces_clears_forces_only() {
        let mut s = two_particle_system();
        {
            let (_, f) = s.split_for_forces();
            f.x[0] = 1.5;
            f.y[1] = -2.5;
            f.z[0] = 3.5;
        }
        assert_eq!(s.forces().x, [1.5, 0.0]);

        s.zero_forces();
        assert_eq!(s.forces().x, [0.0, 0.0]);
        assert_eq!(s.forces().y, [0.0, 0.0]);
        assert_eq!(s.forces().z, [0.0, 0.0]);
        // Positions and velocities untouched.
        assert_eq!(s.positions().x, [1.0, 4.0]);
        assert_eq!(s.velocities().z, [0.3, 0.6]);
    }

    #[test]
    fn split_for_forces_reads_positions_while_writing_forces() {
        let mut s = two_particle_system();
        let sim_box = *s.sim_box();
        let (r, f) = s.split_for_forces();

        // The exact shape of an M1 pair kernel: read r, accumulate into f.
        assert_eq!(r.x.len(), f.x.len());
        let d = sim_box.minimum_image(Vec3::new(r.x[1] - r.x[0], r.y[1] - r.y[0], r.z[1] - r.z[0]));
        f.x[0] += d.x;
        f.x[1] -= d.x;

        assert_relative_eq!(d.x, 3.0);
        assert_relative_eq!(s.forces().x[0], 3.0);
        assert_relative_eq!(s.forces().x[1], -3.0);
    }

    #[test]
    fn box_is_carried_with_the_system() {
        let s = two_particle_system();
        assert_relative_eq!(s.sim_box().volume(), 1000.0);
    }

    #[test]
    #[should_panic(expected = "mass must be finite and positive")]
    fn zero_mass_is_rejected() {
        let mut s = System::with_capacity(1, SimBox::cubic(10.0));
        s.push(Vec3::ZERO, Vec3::ZERO, 0.0, 0.0, 0);
    }
}
