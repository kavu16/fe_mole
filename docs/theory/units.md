# The internal unit system

Implemented in `src/units.rs`. Every constant below is checked against SI by a
test in that file — the tests *derive* each value rather than restating the
literal, so a mistyped constant fails `cargo test` instead of silently biasing
every trajectory the engine ever produces.

## Base units

| Quantity | Unit |
|---|---|
| Length | Å = 10⁻¹⁰ m |
| Mass | amu = 10⁻³ kg·mol⁻¹ / N_A |
| Time | fs = 10⁻¹⁵ s |
| Energy | kcal/mol |
| Temperature | K |
| Charge | e |

Derived: force in kcal/mol/Å, velocity in Å/fs, acceleration in Å/fs².

This combination is the AMBER/CHARMM convention. It is chosen because the
numbers that appear in an MD run land near unity: argon's `ε` is 0.238
kcal/mol, a bond length is a few Å, a timestep is 1–2 fs. Reduced (Lennard-
Jones) units may appear in tests that compare against reduced-unit literature
values, but the conversion happens at the boundary — no core routine sees them.

## SI reference values

The first three are **exact by definition** and will not change:

- `k_B = 1.380649 × 10⁻²³ J/K` (2019 SI redefinition)
- `N_A = 6.02214076 × 10²³ mol⁻¹` (2019 SI redefinition)
- `1 kcal = 4184 J` (definition of the thermochemical calorie)

## Boltzmann constant

Convert per-particle to per-mole, then joules to kilocalories:

```
k_B = 1.380649e-23 J/K · 6.02214076e23 /mol / 4184 J/kcal
    = 8.31446261815324 J/(mol·K) / 4184 J/kcal
    = 1.9872042 × 10⁻³ kcal/(mol·K)
```

`units::BOLTZMANN = 0.0019872041`, which agrees to 8 significant figures.
Tested by `boltzmann_matches_si` at a relative tolerance of `1e-6`.

## Force-to-acceleration conversion

**This is the single most common source of silent error in an MD engine.**
Newton's second law in internal units is *not* `a = f / m`; it needs a
conversion factor, and omitting it produces a simulation that runs, reports
energies, and is wrong by roughly four orders of magnitude.

Take a force of 1 kcal/mol/Å acting on a particle of mass 1 amu.

**Force in newtons, per particle:**

```
1 kcal/mol/Å = 4184 J/mol per Å
             = (4184 / N_A) J per Å           per particle
             = (4184 / N_A) × 10¹⁰ J/m        since 1/Å = 10¹⁰ /m
             = 6.947695 × 10⁻¹¹ N
```

**Mass in kilograms:**

```
1 amu = (10⁻³ kg/mol) / N_A = 1.66053907 × 10⁻²⁷ kg
```

**Acceleration in SI:**

```
a = 6.947695e-11 N / 1.66053907e-27 kg = 4.184 × 10¹⁶ m/s²
```

**Acceleration in internal units:**

```
a = 4.184e16 m/s² × 10¹⁰ Å/m × (10⁻¹⁵ s/fs)²
  = 4.184e16 × 10¹⁰ × 10⁻³⁰
  = 4.184 × 10⁻⁴ Å/fs²
```

So:

```
a [Å/fs²] = 4.184e-4 · f [kcal/mol/Å] / m [amu]
```

### Why the value is exactly 4.184e-4

Not a coincidence, and worth understanding rather than memorising. `N_A`
appears in the numerator (converting kcal/**mol** to per-particle) and in the
denominator (the definition of the amu), so it cancels exactly. What remains
is:

```
4184 × 10¹⁰ × 10¹⁰ × 10⁻³⁰ / 10⁻³ = 4184 × 10⁻⁷ = 4.184 × 10⁻⁴
```

which is exact, because 4184 J/kcal is exact. The test
`force_to_accel_is_exact` therefore asserts a relative tolerance of `1e-12` —
the only error it permits is floating-point round-off in the derivation itself.
That is an unusually tight bound for a physical constant, and it is available
here only because this is a definitional identity rather than a measurement.

### Using it

Prefer `units::acceleration(force, mass)` over multiplying by the constant by
hand. A named function is harder to forget than a multiply, and the failure
mode of forgetting is a plausible-looking wrong answer rather than a crash.

## Sanity checks

Cheap magnitude checks that catch a wrong unit immediately, all in the test
module:

- One argon atom (39.948 amu) under 1 kcal/mol/Å accelerates at
  `1.047 × 10⁻⁵ Å/fs²`.
- At the reference temperature, `k_B·T = 0.1876 kcal/mol` — a bit under
  `ε = k_B × 120 K` for argon, as it must be for a liquid near its triple
  point.
