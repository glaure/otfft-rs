//! The FFT against the direct DFT sum, round trips, scales and refused inputs.

use otfft_rs::{Complex, ComplexFft, Direction, Error, Scale};
use std::f64::consts::PI;

/// The same deterministic signal as the C++ tests of the fork.
fn signal(n: usize) -> Vec<Complex> {
    (0..n)
        .map(|i| Complex::new(((i * 7) % 11) as f64 - 5.0, ((i * 5) % 13) as f64 - 6.0))
        .collect()
}

/// Direct DFT sum of one bin, `sign` -1 forward, +1 inverse, no scale.
fn direct_bin(x: &[Complex], k: usize, sign: f64) -> Complex {
    let n = x.len();
    let mut re = 0.0;
    let mut im = 0.0;
    for (j, v) in x.iter().enumerate() {
        let angle = sign * 2.0 * PI * ((j * k) % n) as f64 / n as f64;
        let w = Complex::expj(angle);
        let p = *v * w;
        re += p.re;
        im += p.im;
    }
    Complex::new(re, im)
}

fn max_diff(a: &[Complex], b: &[Complex]) -> f64 {
    a.iter().zip(b).map(|(p, q)| (*p - *q).abs()).fold(0.0, f64::max)
}

#[test]
fn forward_matches_the_direct_sum() {
    for log_n in 0..=11 {
        let n = 1usize << log_n;
        let input = signal(n);
        let mut x = input.clone();
        ComplexFft::new(n).unwrap().fwd0_alloc(&mut x);
        for (k, got) in x.iter().enumerate() {
            let want = direct_bin(&input, k, -1.0);
            assert!((*got - want).abs() < 1e-9 * n as f64, "n={n} k={k}: {got:?} vs {want:?}");
        }
    }
}

#[test]
fn inverse_matches_the_direct_sum() {
    for log_n in 0..=10 {
        let n = 1usize << log_n;
        let input = signal(n);
        let mut x = input.clone();
        ComplexFft::new(n).unwrap().inv_alloc(&mut x);
        for (k, got) in x.iter().enumerate() {
            let want = direct_bin(&input, k, 1.0);
            assert!((*got - want).abs() < 1e-9 * n as f64, "n={n} k={k}");
        }
    }
}

#[test]
fn fwd_then_inv_returns_the_signal() {
    for log_n in 0..=16 {
        let n = 1usize << log_n;
        let input = signal(n);
        let mut x = input.clone();
        let fft = ComplexFft::new(n).unwrap();
        fft.fwd_alloc(&mut x); // scaled by 1/N
        fft.inv_alloc(&mut x); // not scaled
        assert!(max_diff(&x, &input) < 1e-12 * (n as f64).max(1.0), "n={n}");
    }
}

#[test]
fn the_scales_are_those_of_the_cpp_library() {
    let n = 64;
    let input = signal(n);
    let fft = ComplexFft::new(n).unwrap();

    let mut base = input.clone();
    fft.fwd0_alloc(&mut base);

    for (name, scale, run) in [
        ("fwd", 1.0 / n as f64, ComplexFft::fwd_alloc as fn(&ComplexFft, &mut [Complex])),
        ("fwd0", 1.0, ComplexFft::fwd0_alloc),
        ("fwdu", 1.0 / (n as f64).sqrt(), ComplexFft::fwdu_alloc),
    ] {
        let mut x = input.clone();
        run(&fft, &mut x);
        let want: Vec<Complex> = base.iter().map(|v| Complex::new(v.re * scale, v.im * scale)).collect();
        assert!(max_diff(&x, &want) < 1e-12, "{name}");
    }

    let mut inv_base = input.clone();
    fft.inv_alloc(&mut inv_base);
    for (name, scale, run) in [
        ("inv", 1.0, ComplexFft::inv_alloc as fn(&ComplexFft, &mut [Complex])),
        ("invu", 1.0 / (n as f64).sqrt(), ComplexFft::invu_alloc),
        ("invn", 1.0 / n as f64, ComplexFft::invn_alloc),
    ] {
        let mut x = input.clone();
        run(&fft, &mut x);
        let want: Vec<Complex> = inv_base.iter().map(|v| Complex::new(v.re * scale, v.im * scale)).collect();
        assert!(max_diff(&x, &want) < 1e-12, "{name}");
    }
}

#[test]
fn one_scratch_buffer_serves_many_calls() {
    let n = 256;
    let fft = ComplexFft::new(n).unwrap();
    let mut scratch = vec![Complex::default(); fft.scratch_len()];
    let input = signal(n);
    let mut x = input.clone();
    for _ in 0..3 {
        fft.process(&mut x, &mut scratch, Direction::Forward, Scale::Length).unwrap();
        fft.process(&mut x, &mut scratch, Direction::Inverse, Scale::None).unwrap();
    }
    assert!(max_diff(&x, &input) < 1e-10);
}

#[test]
fn refuses_lengths_and_buffers_it_cannot_handle() {
    assert_eq!(ComplexFft::new(0).unwrap_err(), Error::UnsupportedLength(0));
    assert_eq!(ComplexFft::new(12).unwrap_err(), Error::UnsupportedLength(12));
    let fft = ComplexFft::new(8).unwrap();
    let mut short = vec![Complex::default(); 4];
    let mut scratch = vec![Complex::default(); 8];
    assert_eq!(
        fft.process(&mut short, &mut scratch, Direction::Forward, Scale::None).unwrap_err(),
        Error::WrongBufferLength { expected: 8, got: 4 }
    );
    let mut x = vec![Complex::default(); 8];
    let mut small_scratch = vec![Complex::default(); 2];
    assert_eq!(
        fft.process(&mut x, &mut small_scratch, Direction::Forward, Scale::None).unwrap_err(),
        Error::WrongBufferLength { expected: 8, got: 2 }
    );
}

#[test]
fn a_pure_tone_lands_in_one_bin() {
    let n = 1024;
    let tone: Vec<Complex> = (0..n).map(|j| Complex::expj(2.0 * PI * 50.0 * j as f64 / n as f64)).collect();
    let mut x = tone;
    ComplexFft::new(n).unwrap().fwd_alloc(&mut x);
    for (k, v) in x.iter().enumerate() {
        let want = if k == 50 { 1.0 } else { 0.0 };
        assert!((v.abs() - want).abs() < 1e-9, "bin {k}");
    }
}

#[test]
fn a_plan_can_be_shared_between_threads() {
    fn needs_send_and_sync<T: Send + Sync>() {}
    needs_send_and_sync::<ComplexFft>();

    let fft = std::sync::Arc::new(ComplexFft::new(512).unwrap());
    let input = signal(512);
    let handles: Vec<_> = (0..4)
        .map(|_| {
            let fft = fft.clone();
            let input = input.clone();
            std::thread::spawn(move || {
                let mut x = input.clone();
                fft.fwd_alloc(&mut x);
                fft.inv_alloc(&mut x);
                max_diff(&x, &input)
            })
        })
        .collect();
    for h in handles {
        assert!(h.join().unwrap() < 1e-10);
    }
}

#[test]
fn a_refused_call_leaves_the_data_unchanged() {
    let fft = ComplexFft::new(8).unwrap();
    let input = signal(8);
    let mut x = input.clone();
    let mut small_scratch = vec![Complex::default(); 2];
    assert!(fft.process(&mut x, &mut small_scratch, Direction::Forward, Scale::Length).is_err());
    assert_eq!(x, input);
}

#[test]
fn complex_has_the_layout_of_the_cpp_complex_t() {
    assert_eq!(std::mem::size_of::<Complex>(), 16);
    assert_eq!(std::mem::align_of::<Complex>(), 8);
    let z = Complex::new(1.5, -2.5);
    let raw: [f64; 2] = unsafe { std::mem::transmute(z) };
    assert_eq!(raw, [1.5, -2.5]); // Re first, then Im
}

/// The three remainders of log2(N) modulo 3: radix-8 stages only (2^12), plus one radix-2 stage
/// (2^13), plus one radix-4 stage (2^14).
#[test]
fn large_lengths_of_every_log2_remainder_match_the_direct_sum() {
    for n in [1usize << 12, 1 << 13, 1 << 14] {
        let input = signal(n);
        for (fft_dir, sign) in [(Direction::Forward, -1.0), (Direction::Inverse, 1.0)] {
            let mut x = input.clone();
            let fft = ComplexFft::new(n).unwrap();
            let mut scratch = vec![Complex::default(); fft.scratch_len()];
            fft.process(&mut x, &mut scratch, fft_dir, Scale::None).unwrap();
            // a spread of bins, including the first, the last and the middle ones
            for k in [0, 1, 2, 3, 5, n / 4, n / 2 - 1, n / 2, n / 2 + 1, 3 * n / 4, n - 1] {
                let want = direct_bin(&input, k, sign);
                assert!((x[k] - want).abs() < 1e-9 * n as f64, "n={n} k={k} {fft_dir:?}");
            }
        }
    }
}
