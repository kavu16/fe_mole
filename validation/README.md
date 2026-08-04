# Validation artifacts

Evidence that the acceptance criteria in `MILESTONES.md` actually passed.
These files are committed on purpose — they are the record.

## Rules

1. **Every plot is reproducible from a committed script.** A `.png` with no
   script that regenerates it is not evidence; it is a picture. Commit the
   script alongside it and name them the same.
2. **Every plot states its conditions.** `N`, `dt`, step count, temperature,
   density, cutoff, and the git commit the data came from. A drift curve
   without a `dt` on it means nothing.
3. **Failures stay.** If a run misses its target, the plot stays and `LOG.md`
   records the number. The dead ends are the most useful part of the record
   later, and they are what makes it honest.
4. **Reference data is cited.** Digitised literature curves get a comment
   naming the paper, figure, and how the digitisation was done.

## Expected contents

| Artifact | Milestone | Criterion it supports |
|---|---|---|
| `energy-vs-time.png` | M1 | Relative drift `< 1e-4` over 10⁵ steps |
| `drift-vs-dt.png` | M1 | Drift scales as `O(dt²)` |
| `gr-argon.png` | M2 | `g(r)` overlay against Rahman 1964 |
| `msd-diffusion.png` | M2 | Self-diffusion within ~15% of 2.4e-5 cm²/s |
| `scaling-neighbor-list.png` | M2 | Runtime approximately linear in `N` |
| `kinetic-energy-distribution.png` | M3 | Matches Maxwell–Boltzmann |
| `pme-scaling.png` | M4 | Cost scales as `O(N log N)` |
| `strong-scaling.png`, `weak-scaling.png` | M5 | ≥ 60% efficiency at 8 threads |

Empty at M0 — there is no physics yet, so there is nothing to validate.
