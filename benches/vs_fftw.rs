//! Complex forward transform, `otfft-rs` against RustFFT and FFTW, power-of-two lengths. It
//! follows the method of the C++ benchmark of the fork (`benchmarks/bench.cpp`, `docs/benchmarks.md`
//! there): one forward, unnormalised transform with the plan already built; the in-place
//! transform of `otfft-rs` includes the copy of the input into its work buffer, the out-of-place
//! libraries (RustFFT, FFTW) do not copy; FFTW with `FFTW_MEASURE` and `FFTW_ESTIMATE`, planning
//! limited to 10 s; random input in [-1, 1]; every result is compared with RustFFT before it is
//! timed, and a relative error above 1e-9 stops the run.
//!
//! Needs the system FFTW: `cargo bench --features fftw --bench vs_fftw`. One thread. The
//! harness is criterion, not the batches and medians of the C++ program.

#[path = "support/fftw.rs"]
mod fftw;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use fftw::{FftwPlan, Planning};
use otfft_rs::{Complex, ComplexFft, Direction, Scale};
use rustfft::{FftPlanner, num_complex::Complex64};
use std::hint::black_box;
use std::time::Duration;

/// Random numbers in [-1, 1] (xorshift, the same for every library and every run).
fn random_input(n: usize) -> Vec<(f64, f64)> {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 52) as f64 - 1.0
    };
    (0..n).map(|_| (next(), next())).collect()
}

/// Relative L2 error against the reference; stops the run above 1e-9 like the C++ benchmark.
fn check(name: &str, n: usize, reference: &[Complex64], got: impl Iterator<Item = (f64, f64)>) {
    let (mut diff, mut norm) = (0.0f64, 0.0f64);
    for (r, (re, im)) in reference.iter().zip(got) {
        diff += (r.re - re).powi(2) + (r.im - im).powi(2);
        norm += r.re * r.re + r.im * r.im;
    }
    let error = diff.sqrt() / norm.sqrt();
    assert!(error < 1e-9, "{name} is wrong at n={n}: relative error {error}");
}

fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("vs_fftw");
    group.sample_size(10).warm_up_time(Duration::from_secs(1)).measurement_time(Duration::from_secs(3));
    for log2n in 4..=22 {
        let n = 1usize << log2n;
        let data = random_input(n);
        let rf_input: Vec<Complex64> = data.iter().map(|&(re, im)| Complex64::new(re, im)).collect();

        // RustFFT, out of place; it is also the reference for the other results
        let plan = FftPlanner::<f64>::new().plan_fft_forward(n);
        let mut rf_out = vec![Complex64::default(); n];
        let mut rf_scratch = vec![Complex64::default(); plan.get_outofplace_scratch_len()];
        let mut rf_in = rf_input.clone();
        plan.process_outofplace_with_scratch(&mut rf_in, &mut rf_out, &mut rf_scratch);
        let reference = rf_out.clone();

        // otfft-rs, in place: copy the input, then transform
        let fft = ComplexFft::new(n).unwrap();
        let source: Vec<Complex> = data.iter().map(|&(re, im)| Complex::new(re, im)).collect();
        let mut work = source.clone();
        let mut scratch = vec![Complex::default(); fft.scratch_len()];
        fft.process(&mut work, &mut scratch, Direction::Forward, Scale::None).unwrap();
        check("otfft-rs", n, &reference, work.iter().map(|v| (v.re, v.im)));

        let mut measure = FftwPlan::new(n, Planning::Measure);
        measure.set_input(data.iter().copied());
        measure.execute();
        check("fftw-measure", n, &reference, measure.output());

        let mut estimate = FftwPlan::new(n, Planning::Estimate);
        estimate.set_input(data.iter().copied());
        estimate.execute();
        check("fftw-estimate", n, &reference, estimate.output());

        group.bench_with_input(BenchmarkId::new("otfft-rs", n), &n, |b, _| {
            b.iter(|| {
                work.copy_from_slice(black_box(&source));
                fft.process(&mut work, &mut scratch, Direction::Forward, Scale::None).unwrap()
            })
        });
        group.bench_with_input(BenchmarkId::new("rustfft", n), &n, |b, _| {
            b.iter(|| {
                plan.process_outofplace_with_scratch(
                    black_box(&mut rf_in),
                    &mut rf_out,
                    &mut rf_scratch,
                )
            })
        });
        group.bench_with_input(BenchmarkId::new("fftw-measure", n), &n, |b, _| {
            b.iter(|| measure.execute())
        });
        group.bench_with_input(BenchmarkId::new("fftw-estimate", n), &n, |b, _| {
            b.iter(|| estimate.execute())
        });
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
