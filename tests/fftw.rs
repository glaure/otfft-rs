//! The FFTW binding of the benchmark (`benches/support/fftw.rs`) and otfft-rs against FFTW.
//! Runs only with `--features fftw` and a system `libfftw3`. The binding has `unsafe` code; this
//! test runs it (AGENTS.md rule 4).
#![cfg(feature = "fftw")]

#[path = "../benches/support/fftw.rs"]
mod fftw;

use fftw::{FftwPlan, Planning};
use otfft_rs::{Complex, ComplexFft};
use std::f64::consts::PI;

fn input(n: usize) -> Vec<(f64, f64)> {
    (0..n).map(|i| (((i * 7) % 11) as f64 - 5.0, ((i * 5) % 13) as f64 - 6.0)).collect()
}

#[test]
fn the_binding_computes_the_direct_sum() {
    // Checks sign, layout and scale of the binding itself, so that "equal to FFTW" means
    // something in the next test.
    let n = 64;
    let data = input(n);
    for planning in [Planning::Estimate, Planning::Measure] {
        let mut plan = FftwPlan::new(n, planning);
        plan.set_input(data.iter().copied());
        plan.execute();
        for (k, (re, im)) in plan.output().enumerate() {
            let (mut wr, mut wi) = (0.0, 0.0);
            for (j, &(xr, xi)) in data.iter().enumerate() {
                let a = -2.0 * PI * ((j * k) % n) as f64 / n as f64;
                wr += xr * a.cos() - xi * a.sin();
                wi += xr * a.sin() + xi * a.cos();
            }
            assert!((re - wr).abs() < 1e-9 && (im - wi).abs() < 1e-9, "{planning:?} k={k}");
        }
    }
}

#[test]
fn otfft_rs_equals_fftw() {
    for log_n in 0..=16 {
        let n = 1usize << log_n;
        let data = input(n);
        let mut plan = FftwPlan::new(n, Planning::Estimate);
        plan.set_input(data.iter().copied());
        plan.execute();

        let mut x: Vec<Complex> = data.iter().map(|&(re, im)| Complex::new(re, im)).collect();
        ComplexFft::new(n).unwrap().fwd0_alloc(&mut x);

        let (mut diff, mut norm) = (0.0f64, 0.0f64);
        for (got, (re, im)) in x.iter().zip(plan.output()) {
            diff += (got.re - re).powi(2) + (got.im - im).powi(2);
            norm += re * re + im * im;
        }
        assert!(diff.sqrt() <= 1e-12 * norm.sqrt().max(1.0), "n={n}: {diff} vs {norm}");
    }
}
