//! Complex forward transform, `otfft-rs` against RustFFT, power-of-two lengths.
//!
//! Both transforms work in place with a caller-owned scratch buffer and no scaling, so the
//! timings compare the transform and not the allocation. `realfft` joins in with the real FFT
//! (README step 5).

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use otfft_rs::{ComplexFft, Direction, Scale};
use rustfft::{FftPlanner, num_complex::Complex64};
use std::hint::black_box;

fn input(n: usize) -> Vec<(f64, f64)> {
    (0..n).map(|i| ((i as f64 * 0.37).sin(), (i as f64 * 0.11).cos())).collect()
}

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("complex_fwd");
    for log2n in [4, 6, 8, 10, 11, 12, 14, 16, 17, 18, 20] {
        let n = 1usize << log2n;
        group.throughput(Throughput::Elements(n as u64));
        let data = input(n);

        let fft = ComplexFft::new(n).unwrap();
        let mut x: Vec<_> = data.iter().map(|&(re, im)| otfft_rs::Complex::new(re, im)).collect();
        let mut scratch = vec![otfft_rs::Complex::default(); fft.scratch_len()];
        group.bench_with_input(BenchmarkId::new("otfft-rs", n), &n, |b, _| {
            b.iter(|| {
                fft.process(black_box(&mut x), &mut scratch, Direction::Forward, Scale::None)
                    .unwrap()
            })
        });

        let plan = FftPlanner::<f64>::new().plan_fft_forward(n);
        let mut y: Vec<_> = data.iter().map(|&(re, im)| Complex64::new(re, im)).collect();
        let mut scratch = vec![Complex64::default(); plan.get_inplace_scratch_len()];
        group.bench_with_input(BenchmarkId::new("rustfft", n), &n, |b, _| {
            b.iter(|| plan.process_with_scratch(black_box(&mut y), &mut scratch))
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
