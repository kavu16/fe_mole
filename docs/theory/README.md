# Theory notes

Distilled notes on the physics and numerics behind the engine, using the
**same notation as the code**. These exist so that a derivation and its
implementation can be read side by side; they are not a substitute for the
cited sources, and they deliberately do not restate textbook derivations that
add nothing.

## Index

| Note | Covers | Milestone |
|---|---|---|
| [`units.md`](units.md) | Internal unit system; derivation of `BOLTZMANN` and `FORCE_TO_ACCEL` | M0 |
| _`lennard-jones.md`_ | LJ potential, cutoff treatment, shift vs. switch | M1 |
| _`integrators.md`_ | Velocity Verlet as a symplectic splitting; why drift is `O(dt²)` | M1 |
| _`neighbor-lists.md`_ | Cell lists, Verlet lists, skin distance, rebuild criterion | M2 |
| _`thermostats.md`_ | Langevin dynamics; BAOAB vs. naive splitting | M3 |
| _`pme.md`_ | Ewald splitting; smooth PME with B-splines | M4 |

Italicised notes are not written yet.

## Notation

Shared by the notes and the code. Where a symbol appears in an identifier, the
identifier uses the ASCII form in the third column.

| Symbol | Meaning | In code | Unit |
|---|---|---|---|
| `N` | number of particles | `n`, `System::len()` | — |
| `rᵢ` | position of particle `i` (**unwrapped**) | `rx[i]`, `ry[i]`, `rz[i]` | Å |
| `vᵢ` | velocity of particle `i` | `vx[i]`, `vy[i]`, `vz[i]` | Å/fs |
| `fᵢ` | force on particle `i` | `fx[i]`, `fy[i]`, `fz[i]` | kcal/mol/Å |
| `rᵢⱼ` | minimum-image displacement `rⱼ − rᵢ` | `d` | Å |
| `mᵢ` | mass of particle `i` | `mass[i]` | amu |
| `qᵢ` | charge of particle `i` | `charge[i]` | e |
| `L` | box side lengths | `SimBox::lengths()` | Å |
| `V` | box volume | `SimBox::volume()` | Å³ |
| `ρ` | mass density | `density` | g/cm³ |
| `σ` | LJ length parameter | `sigma` | Å |
| `ε` | LJ energy parameter | `epsilon` | kcal/mol |
| `r_c` | cutoff radius | `r_cut` | Å |
| `r_s` | neighbour-list skin | `r_skin` | Å |
| `dt` | integration timestep | `dt` | fs |
| `T` | temperature | `temperature` | K |
| `k_B` | Boltzmann constant | `units::BOLTZMANN` | kcal/(mol·K) |
| `U` | potential energy (one pair, or summed) | `energy` | kcal/mol |
| `E` | total energy, `KE + U` | `energy_total` | kcal/mol |

Two conventions are worth stating because they are easy to get backwards:

- `rᵢⱼ = rⱼ − rᵢ`, so the force on `i` from `j` points along `+rᵢⱼ` when the
  interaction is attractive.
- Positions are **never wrapped in place**. `rᵢ` accumulates without bound;
  wrapping happens only at output. See `system::System` for why.

## Reference system

Everything in the milestone list is validated against Rahman's 1964 liquid
argon simulation:

- `N = 864`, `T = 94.4 K`, `ρ = 1.374 g/cm³`
- `σ = 3.4 Å`, `ε/k_B = 120 K`
- Cubic box, `L ≈ 34.68 Å` (computed by `SimBox::from_density`)

A. Rahman, "Correlations in the Motion of Atoms in Liquid Argon",
*Phys. Rev.* **136**, A405 (1964).
