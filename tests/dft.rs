//! The FFT against the direct DFT sum, round trips, scales and refused inputs.

use otfft_rs::{Isa, Complex, ComplexFft, Direction, Error, Kernel, Scale};
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

/// The remainders of log2(N) for every class: modulo 4 for the radix-16 class (radix-16 stages
/// only at 2^12, plus one radix-2 stage at 2^13, one radix-4 stage at 2^14, one radix-8 stage at
/// 2^15), modulo 3 for the radix-8 class, both parities for the radix-4 class. Both directions.
#[test]
fn large_lengths_of_every_log2_remainder_match_the_direct_sum() {
    for &kernel in Kernel::ALL {
        for n in [1usize << 12, 1 << 13, 1 << 14, 1 << 15] {
            let input = signal(n);
            for (fft_dir, sign) in [(Direction::Forward, -1.0), (Direction::Inverse, 1.0)] {
                let mut x = input.clone();
                let fft = ComplexFft::with_kernel(n, kernel).unwrap();
                let mut scratch = vec![Complex::default(); fft.scratch_len()];
                fft.process(&mut x, &mut scratch, fft_dir, Scale::None).unwrap();
                // a spread of bins, including the first, the last and the middle ones
                for k in [0, 1, 2, 3, 5, n / 4, n / 2 - 1, n / 2, n / 2 + 1, 3 * n / 4, n - 1] {
                    let want = direct_bin(&input, k, sign);
                    assert!(
                        (x[k] - want).abs() < 1e-9 * n as f64,
                        "{kernel:?} n={n} k={k} {fft_dir:?}"
                    );
                }
            }
        }
    }
}

/// Every algorithm class matches the direct sum in every bin, both directions, for the lengths 1
/// to 2^11 (all the last-stage cases of every class).
#[test]
fn every_kernel_matches_the_direct_sum() {
    for &kernel in Kernel::ALL {
        for log_n in 0..=11 {
            let n = 1usize << log_n;
            let input = signal(n);
            for (dir, sign) in [(Direction::Forward, -1.0), (Direction::Inverse, 1.0)] {
                let fft = ComplexFft::with_kernel(n, kernel).unwrap();
                assert_eq!(fft.kernel(), kernel);
                let mut x = input.clone();
                let mut scratch = vec![Complex::default(); fft.scratch_len()];
                fft.process(&mut x, &mut scratch, dir, Scale::None).unwrap();
                for (k, got) in x.iter().enumerate() {
                    let want = direct_bin(&input, k, sign);
                    assert!(
                        (*got - want).abs() < 1e-9 * n as f64,
                        "{kernel:?} n={n} k={k} {dir:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn new_uses_the_default_kernel_and_with_kernel_refuses_bad_lengths() {
    for log_n in 0..=12 {
        let n = 1usize << log_n;
        assert_eq!(ComplexFft::new(n).unwrap().kernel(), Kernel::default_for(n));
    }
    for &kernel in Kernel::ALL {
        for n in [0, 3, 12, 1000] {
            assert_eq!(ComplexFft::with_kernel(n, kernel).unwrap_err(), Error::UnsupportedLength(n));
        }
    }
}

/// Every instruction set that this CPU has gives the result of the direct sum in every bin, and
/// the same bits as the scalar stages, for every class, both directions, lengths 1 to 2^11 (the
/// AVX2 stages need a stride of 2 or more, so the short lengths cover the scalar fallback inside
/// an AVX2 plan).
#[test]
fn every_isa_matches_the_direct_sum_and_the_scalar_bits() {
    for &kernel in Kernel::ALL {
        for &isa in Isa::ALL.iter().filter(|i| i.is_available()) {
            for log_n in 0..=11 {
                let n = 1usize << log_n;
                let input = signal(n);
                for (dir, sign) in [(Direction::Forward, -1.0), (Direction::Inverse, 1.0)] {
                    let run = |isa| {
                        let fft = ComplexFft::with_kernel_isa(n, kernel, isa).unwrap();
                        assert_eq!(fft.isa(), isa);
                        let mut x = input.clone();
                        let mut scratch = vec![Complex::default(); fft.scratch_len()];
                        fft.process(&mut x, &mut scratch, dir, Scale::None).unwrap();
                        x
                    };
                    let (x, scalar) = (run(isa), run(Isa::Scalar));
                    assert_eq!(x, scalar, "{kernel:?} {isa:?} n={n} {dir:?}: not the scalar bits");
                    for (k, got) in x.iter().enumerate() {
                        let want = direct_bin(&input, k, sign);
                        assert!(
                            (*got - want).abs() < 1e-9 * n as f64,
                            "{kernel:?} {isa:?} n={n} k={k} {dir:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn an_unavailable_isa_is_refused_and_detect_is_available() {
    assert!(Isa::Scalar.is_available());
    assert!(Isa::detect().is_available());
    for &isa in Isa::ALL {
        let plan = ComplexFft::with_kernel_isa(64, Kernel::Dif4, isa);
        if isa.is_available() {
            assert_eq!(plan.unwrap().isa(), isa);
        } else {
            assert_eq!(plan.unwrap_err(), Error::IsaUnavailable(isa));
        }
    }
    assert_eq!(ComplexFft::new(64).unwrap().isa(), Isa::detect());
}
