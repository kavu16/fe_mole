//! The internal unit system, as named constants.
//!
//! | Quantity | Unit |
//! |---|---|
//! | Length | Å |
//! | Mass | amu |
//! | Time | fs |
//! | Energy | kcal/mol |
//! | Temperature | K |
//! | Charge | e |
//!
//! Derived: force is kcal/mol/Å, velocity is Å/fs, acceleration is Å/fs².
//!
//! Every constant here is checked against SI in the tests below — the tests
//! *derive* the value rather than restating the literal, so a typo in a
//! constant fails the build rather than silently biasing every trajectory.
//! The derivations are written out in `docs/theory/units.md`.

/// Avogadro constant, mol⁻¹. Exact by the 2019 SI redefinition.
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Boltzmann constant in internal units, kcal/(mol·K).
///
/// Equal to `k_B[SI] · N_A / 4184`.
pub const BOLTZMANN: f64 = 0.0019872041;

/// Converts a force over a mass into an acceleration:
/// `a [Å/fs²] = FORCE_TO_ACCEL · f [kcal/mol/Å] / m [amu]`.
///
/// **This factor is not optional.** Omitting it produces a simulation that
/// looks entirely plausible — particles move, energy is reported — and is
/// wrong by four orders of magnitude. Prefer [`acceleration`] over using this
/// constant directly; a named function is harder to forget than a multiply.
///
/// The value is *exact*: the Avogadro factor cancels between "per mole" in the
/// energy unit and the definition of the amu, and the thermochemical calorie
/// is defined as exactly 4184 J.
pub const FORCE_TO_ACCEL: f64 = 4.184e-4;

/// Converts `m · v²` from internal units (amu·Å²/fs²) to energy in kcal/mol:
/// `E = MASS_VELOCITY_SQ_TO_ENERGY · m · v²`.
///
/// Needed by every kinetic-energy and temperature calculation. Exactly
/// `1 / FORCE_TO_ACCEL`, and necessarily so: the factor that turns a force
/// into an acceleration has to turn the resulting work back into an energy,
/// or kinetic and potential energy land on different scales and their sum is
/// meaningless.
///
/// Omitting it inflates kinetic energy by a factor of ~2390, which surfaces as
/// an absurd temperature rather than a slow drift — one of the few unit errors
/// in MD that fails loudly.
pub const MASS_VELOCITY_SQ_TO_ENERGY: f64 = 1.0 / FORCE_TO_ACCEL;

/// Acceleration [Å/fs²] of a particle of mass `mass` [amu] under a force
/// component `force` [kcal/mol/Å].
///
/// This is the only sanctioned way to turn a force into an acceleration.
#[inline]
#[must_use]
pub fn acceleration(force: f64, mass: f64) -> f64 {
    FORCE_TO_ACCEL * force / mass
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // SI reference values. The first three are exact by definition since the
    // 2019 SI redefinition and the definition of the thermochemical calorie;
    // they are not measured quantities and will not change.
    const K_B_SI: f64 = 1.380_649e-23; // J/K
    const JOULES_PER_KCAL: f64 = 4184.0; // exact, thermochemical
    const ANGSTROM_PER_METRE: f64 = 1e10;
    const FS_PER_SECOND: f64 = 1e15;
    /// 1 amu = (molar mass constant) / N_A = 10⁻³ kg/mol / N_A.
    const KG_PER_AMU: f64 = 1e-3 / AVOGADRO;

    #[test]
    fn boltzmann_matches_si() {
        // k_B [J/K] -> [J/(mol·K)] -> [kcal/(mol·K)]
        let derived = K_B_SI * AVOGADRO / JOULES_PER_KCAL;
        assert_relative_eq!(BOLTZMANN, derived, max_relative = 1e-6);
    }

    #[test]
    fn force_to_accel_is_exact() {
        // A force of 1 kcal/mol/Å, in newtons on a single particle:
        //   (4184 J/mol) / N_A per Å, and 1/Å = 1e10 /m.
        let force_si = JOULES_PER_KCAL / AVOGADRO * ANGSTROM_PER_METRE;
        let accel_si = force_si / KG_PER_AMU; // m/s²
        let derived = accel_si * ANGSTROM_PER_METRE / (FS_PER_SECOND * FS_PER_SECOND);

        // This is an analytical identity, not a measurement: the only error
        // permitted is floating-point round-off in the derivation above.
        assert_relative_eq!(FORCE_TO_ACCEL, derived, max_relative = 1e-12);
    }

    #[test]
    fn mass_velocity_sq_to_energy_matches_si() {
        // One amu moving at 1 Å/fs, taken to J and then to kcal/mol.
        let v_si = FS_PER_SECOND / ANGSTROM_PER_METRE; // 1 Å/fs in m/s
        let energy_si = KG_PER_AMU * v_si * v_si;
        let derived = energy_si * AVOGADRO / JOULES_PER_KCAL;
        assert_relative_eq!(MASS_VELOCITY_SQ_TO_ENERGY, derived, max_relative = 1e-12);
    }

    #[test]
    fn energy_and_acceleration_conversions_are_reciprocal() {
        // Not a coincidence: work done by a force must come back out as an
        // energy on the same scale the force was defined on.
        assert_relative_eq!(
            MASS_VELOCITY_SQ_TO_ENERGY * FORCE_TO_ACCEL,
            1.0,
            max_relative = 1e-15
        );
    }

    #[test]
    fn thermal_speed_of_argon() {
        // Magnitude check against a hand calculation: v_rms = sqrt(3kT/m) for
        // argon at the reference temperature is 242.8 m/s, i.e. 2.428e-3 Å/fs.
        let v_rms_sq = 3.0 * BOLTZMANN * 94.4 / MASS_VELOCITY_SQ_TO_ENERGY / 39.948;
        assert_relative_eq!(v_rms_sq.sqrt(), 2.428e-3, max_relative = 1e-3);
    }

    #[test]
    fn accelerating_argon() {
        // 1 kcal/mol/Å on one argon atom (39.948 amu). By hand:
        // 4.184e-4 / 39.948 = 1.0474e-5 Å/fs².
        let a = acceleration(1.0, 39.948);
        assert_relative_eq!(a, 1.047_361e-5, max_relative = 1e-6);
    }

    #[test]
    fn acceleration_is_linear_in_force_and_inverse_in_mass() {
        assert_relative_eq!(acceleration(2.0, 10.0), 2.0 * acceleration(1.0, 10.0));
        assert_relative_eq!(acceleration(1.0, 20.0), 0.5 * acceleration(1.0, 10.0));
    }

    #[test]
    fn thermal_energy_at_reference_temperature() {
        // Rahman 1964 liquid argon runs at T = 94.4 K; ε/k_B = 120 K.
        // k_B·T should therefore be a bit under ε — a cheap check that
        // BOLTZMANN is in kcal/mol/K and not J/K or reduced units.
        let kt = BOLTZMANN * 94.4;
        let epsilon = BOLTZMANN * 120.0;
        assert!(kt < epsilon);
        assert_relative_eq!(kt / epsilon, 94.4 / 120.0, max_relative = 1e-12);
        // ~0.1876 kcal/mol; if this were off by 10³ the unit is wrong.
        assert!((0.18..0.19).contains(&kt), "k_B·T = {kt} kcal/mol");
    }
}
