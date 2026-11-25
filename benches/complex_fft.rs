//! Complex forward transform, every algorithm class of `otfft-rs` against RustFFT, power-of-two
//! lengths.
//!
//! Both transforms work in place with a caller-owned scratch buffer and no scaling, so the
//! timings compare the transform and not the allocation. `realfft` joins in with the real FFT
//! (README step 5).

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use otfft_rs::{ComplexFft, Direction, Isa, Kernel, Scale};
use rustfft::{FftPlanner, num_complex::Complex64};
use std::hint::black_box;

fn input(n: usize) -> Vec<(f64, f64)> {
    (0..n).map(|i| ((i as f64 * 0.37).sin(), (i as f64 * 0.11).cos())).collect()
}

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("complex_fwd");
    for log2n in [4, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20] {
        let n = 1usize << log2n;
        group.throughput(Throughput::Elements(n as u64));
        let data = input(n);

        // Scalar stages for every class, so the numbers of B1 to B9 stay comparable, and the AVX2
        // stages (named `-Avx2`) for the classes that have them.
        for &kernel in Kernel::ALL {
            for &isa in Isa::ALL {
                if !isa.is_available() || (isa != Isa::Scalar && !kernel.has_avx2()) {
                    continue;
                }
                let fft = ComplexFft::with_kernel_isa(n, kernel, isa).unwrap();
                let mut x: Vec<_> =
                    data.iter().map(|&(re, im)| otfft_rs::Complex::new(re, im)).collect();
                let mut scratch = vec![otfft_rs::Complex::default(); fft.scratch_len()];
                let name = if isa == Isa::Scalar {
                    format!("otfft-rs-{kernel:?}")
                } else {
                    format!("otfft-rs-{kernel:?}-{isa:?}")
                };
                group.bench_with_input(BenchmarkId::new(name, n), &n, |b, _| {
                    b.iter(|| {
                        fft.process(black_box(&mut x), &mut scratch, Direction::Forward, Scale::None)
                            .unwrap()
                    })
                });
            }
        }

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
