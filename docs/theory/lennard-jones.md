# The Lennard-Jones potential

*Skeleton — fill in as checkpoint 2 is implemented. Use the notation in
[`README.md`](README.md); the code in `src/forces/lj.rs` should read as a
direct transcription of what ends up here.*

## The potential

State the 12-6 form and what each term is doing physically — which one is the
Pauli repulsion, which is the induced-dipole attraction, and which of the two
has a defensible physical origin. (Only one does. The other is chosen for a
reason worth writing down.)

Note the two equivalent parameterisations (`σ`/`ε` versus `r_min`/`ε`) and
which one this code uses, since mixing them up is a factor of `2^(1/6)`.

## The force

Derive `f(r) = −dV/dr`. Keep the algebra, not just the result — the exponents
are where sign and factor errors live.

Then note what the code actually computes and why: the pair loop wants
`−V'(r)/r` so it can multiply by the displacement vector directly, avoiding a
square root per pair. Record what that costs in clarity and what it buys.

## Landmarks

Three values worth deriving, because they are the cheapest possible checks on
an implementation:

- where `V(r) = 0`
- where `V(r)` is minimised, and its value there
- where `f(r)` is maximally attractive

## Cutoff treatment

Why a cutoff is needed at all (cost, and the minimum image constraint
`r_c ≤ L/2`).

Then the central point: an unmodified cutoff leaves `V` discontinuous at
`r_c`, and every pair crossing it injects energy — this shows up as monotonic
heating in an NVE run, not as random drift. Explain the shifted form
`V_shift(r) = V(r) − V(r_c)`, and answer:

- What does shifting do to the force? (Work out `d/dr` of a constant. This is
  the point of the whole exercise.)
- Is the *force* continuous at `r_c` under this scheme? If not, what is the
  remaining error, and what would fix it?
- What is thrown away by truncating, and what is the standard correction for
  the discarded tail in energy and pressure?

Record which variant this engine implements and why, per `CLAUDE.md`'s rule on
naming the variant.

## Reduced units

The literature quotes liquid-argon results in reduced units (`r* = r/σ`,
`T* = k_B T/ε`, `ρ* = ρσ³`). Give the conversions, since M2 compares `g(r)`
against Rahman 1964. Conversion happens at the test boundary — the engine core
stays in Å / kcal/mol / K.

## Parameters

| Species | `σ` (Å) | `ε/k_B` (K) | `ε` (kcal/mol) | Source |
|---|---|---|---|---|
| Ar | 3.4 | 120 | | Rahman 1964 |

## Sources

- J. E. Lennard-Jones, *Proc. Phys. Soc.* **43**, 461 (1931)
- A. Rahman, *Phys. Rev.* **136**, A405 (1964)
- Allen & Tildesley, *Computer Simulation of Liquids*, ch. 1 and §2.8
  (cutoffs, shifted potentials, tail corrections)
