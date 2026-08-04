//! Benchmark for [`SimBox::minimum_image`].
//!
//! M0 asks for one trivial benchmark, purely to prove the harness produces a
//! number. This one was chosen because it will stay useful: the minimum image
//! convention is applied to every pair displacement in every force
//! evaluation, so this is the innermost operation of M1's O(N²) loop, M2's
//! neighbour lists and M4's real-space Ewald sum. Whatever the per-pair cost
//! turns out to be, it is a floor on all of them.
//!
//! Two variants are timed so the numbers are comparable later:
//!
//! * `round` — the form the engine actually uses. Correct for a displacement
//!   of any magnitude.
//! * `branch` — the cheaper conditional form, valid only for `|d| < 3L/2`.
//!   Benchmarked, **not** exported, and not used anywhere in the engine. It is
//!   here to quantify what the correct version costs, so that any future
//!   argument to switch has a number attached to it.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fe_mole::geometry::{SimBox, Vec3};
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;

/// Number of displacements per timed iteration.
const N_PAIRS: usize = 10_000;

/// Fixed seed: benchmark inputs must be identical run to run, or a "5% faster"
/// result is indistinguishable from different random data.
const SEED: u64 = 0xB0DE_1234;

/// Generates displacements spanning several box widths, struct-of-arrays, to
/// match how the force loop will actually feed this routine.
fn displacements(sim_box: &SimBox, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let l = sim_box.lengths();
    let mut rng = Pcg64Mcg::seed_from_u64(SEED);
    let mut dx = Vec::with_capacity(n);
    let mut dy = Vec::with_capacity(n);
    let mut dz = Vec::with_capacity(n);
    for _ in 0..n {
        dx.push(rng.random_range(-1.4 * l.x..1.4 * l.x));
        dy.push(rng.random_range(-1.4 * l.y..1.4 * l.y));
        dz.push(rng.random_range(-1.4 * l.z..1.4 * l.z));
    }
    (dx, dy, dz)
}

/// The branch form, for comparison only. Valid only while `|d| < 3L/2`, which
/// is why the input range above stays inside ±1.4 L.
#[inline]
fn minimum_image_branch(d: f64, l: f64) -> f64 {
    if d > 0.5 * l {
        d - l
    } else if d < -0.5 * l {
        d + l
    } else {
        d
    }
}

fn bench_minimum_image(c: &mut Criterion) {
    // Rahman 1964 reference box: 864 argon atoms at 1.374 g/cm³.
    let sim_box = SimBox::from_density(864, 39.948, 1.374);
    let l = sim_box.lengths();
    let (dx, dy, dz) = displacements(&sim_box, N_PAIRS);

    let mut group = c.benchmark_group("minimum_image");
    group.throughput(criterion::Throughput::Elements(N_PAIRS as u64));

    group.bench_function("round", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for i in 0..N_PAIRS {
                let d = sim_box.minimum_image(Vec3::new(dx[i], dy[i], dz[i]));
                acc += d.norm_squared();
            }
            black_box(acc)
        });
    });

    group.bench_function("branch", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for i in 0..N_PAIRS {
                let x = minimum_image_branch(dx[i], l.x);
                let y = minimum_image_branch(dy[i], l.y);
                let z = minimum_image_branch(dz[i], l.z);
                acc += x * x + y * y + z * z;
            }
            black_box(acc)
        });
    });

    group.finish();
}

criterion_group!(benches, bench_minimum_image);
criterion_main!(benches);
