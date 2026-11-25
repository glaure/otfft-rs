//! AVX2 stages (README step 4, D17). One `__m256d` holds two `Complex` values.
//!
//! Every function here is `unsafe` for one reason: it has the `avx2` target feature, so the caller
//! must have checked `is_x86_feature_detected!("avx2")`. `ComplexFft` does that once, in
//! `Isa::is_available`, and only a plan with `Isa::Avx2` calls into this module. The memory
//! accesses are inside the slices: each function asserts the length of its buffers on entry, and
//! the index arithmetic below stays under that length (comment at each access).
//!
//! The arithmetic is the same additions and multiplications as the scalar stages, in the same
//! order and without fused multiply-add, so the result equals the scalar one bit for bit.

use crate::Complex;
use std::arch::x86_64::*;

/// `a * w` for two complex values in `a` and the factor `w` (real and imaginary part broadcast
/// in `wr` and `wi`).
#[inline]
#[target_feature(enable = "avx2")]
fn mul(a: __m256d, wr: __m256d, wi: __m256d) -> __m256d {
    let swapped = _mm256_permute_pd::<0b0101>(a); // [im, re] of each value
    // even lanes: re * wr - im * wi, odd lanes: im * wr + re * wi
    _mm256_addsub_pd(_mm256_mul_pd(a, wr), _mm256_mul_pd(swapped, wi))
}

/// Multiplication by -i (forward) or +i (inverse): `(re, im)` to `(im, -re)` or `(-im, re)`.
#[inline]
#[target_feature(enable = "avx2")]
fn mul_j(z: __m256d, inverse: bool) -> __m256d {
    let swapped = _mm256_permute_pd::<0b0101>(z);
    let sign = if inverse { _mm256_set_pd(0.0, -0.0, 0.0, -0.0) } else { _mm256_set_pd(-0.0, 0.0, -0.0, 0.0) };
    _mm256_xor_pd(swapped, sign)
}

#[inline]
#[target_feature(enable = "avx2")]
fn load(p: *const Complex) -> __m256d {
    // SAFETY: the callers pass a pointer to two `Complex` (32 bytes) inside a slice.
    unsafe { _mm256_loadu_pd(p as *const f64) }
}

#[inline]
#[target_feature(enable = "avx2")]
fn store(p: *mut Complex, v: __m256d) {
    // SAFETY: the callers pass a pointer to two `Complex` (32 bytes) inside a slice.
    unsafe { _mm256_storeu_pd(p as *mut f64, v) }
}

/// The 4-point DFT of `a, b, c, d` without twiddles, in the order of the scalar `radix4_stage`.
#[inline]
#[target_feature(enable = "avx2")]
fn dft4(a: __m256d, b: __m256d, c: __m256d, d: __m256d, inverse: bool) -> [__m256d; 4] {
    let apc = _mm256_add_pd(a, c);
    let amc = _mm256_sub_pd(a, c);
    let bpd = _mm256_add_pd(b, d);
    let jbmd = mul_j(_mm256_sub_pd(b, d), inverse);
    [
        _mm256_add_pd(apc, bpd),
        _mm256_add_pd(amc, jbmd),
        _mm256_sub_pd(apc, bpd),
        _mm256_sub_pd(amc, jbmd),
    ]
}

/// The AVX2 form of `ComplexFft::radix4_stage`: current length `n` (at least 8), stride `s`, with
/// `s` a power of two and at least 2. `twiddle` is the table of the plan (`total` entries).
///
/// # Safety
/// The CPU must support AVX2.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn radix4_stage(
    src: &[Complex],
    dst: &mut [Complex],
    twiddle: &[Complex],
    n: usize,
    s: usize,
    inverse: bool,
) {
    let total = n * s; // the length of the plan
    assert!(s >= 2 && s.is_multiple_of(2) && n >= 8 && twiddle.len() == total);
    assert!(src.len() >= total && dst.len() >= total);
    let m = n / 4;
    let (sp, dp) = (src.as_ptr(), dst.as_mut_ptr());
    let plan_step = total / n; // twiddle index of `exp(-2 pi i / n)`
    for p in 0..m {
        let tw = |k: usize| {
            let w = twiddle[k * p * plan_step];
            if inverse { w.conj() } else { w }
        };
        let (w1, w2, w3) = (tw(1), tw(2), tw(3));
        let (w1r, w1i) = (_mm256_set1_pd(w1.re), _mm256_set1_pd(w1.im));
        let (w2r, w2i) = (_mm256_set1_pd(w2.re), _mm256_set1_pd(w2.im));
        let (w3r, w3i) = (_mm256_set1_pd(w3.re), _mm256_set1_pd(w3.im));
        for q in (0..s).step_by(2) {
            // Indices are below `s * n = total`: `q + 1 < s`, `p + 3 * m < n`.
            // SAFETY: in bounds by the assertion on the lengths, see above.
            unsafe {
                let a = load(sp.add(q + s * p));
                let b = load(sp.add(q + s * (p + m)));
                let c = load(sp.add(q + s * (p + 2 * m)));
                let d = load(sp.add(q + s * (p + 3 * m)));
                let [o0, o1, o2, o3] = dft4(a, b, c, d, inverse);
                store(dp.add(q + s * (4 * p)), o0);
                store(dp.add(q + s * (4 * p + 1)), mul(o1, w1r, w1i));
                store(dp.add(q + s * (4 * p + 2)), mul(o2, w2r, w2i));
                store(dp.add(q + s * (4 * p + 3)), mul(o3, w3r, w3i));
            }
        }
    }
}

/// The AVX2 form of `ComplexFft::radix4_last`: `s` a power of two and at least 2.
///
/// # Safety
/// The CPU must support AVX2.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn radix4_last(src: &[Complex], dst: &mut [Complex], s: usize, inverse: bool) {
    assert!(s >= 2 && s.is_multiple_of(2) && src.len() >= 4 * s && dst.len() >= 4 * s);
    let (sp, dp) = (src.as_ptr(), dst.as_mut_ptr());
    for q in (0..s).step_by(2) {
        // SAFETY: `q + 1 < s`, so every index is below `4 * s`, the asserted length.
        unsafe {
            let o = dft4(
                load(sp.add(q)),
                load(sp.add(q + s)),
                load(sp.add(q + 2 * s)),
                load(sp.add(q + 3 * s)),
                inverse,
            );
            for (k, v) in o.into_iter().enumerate() {
                store(dp.add(q + k * s), v);
            }
        }
    }
}

/// Multiplication of both values of `z` by the real `k`, as `Complex::scale`.
#[inline]
#[target_feature(enable = "avx2")]
fn scale(z: __m256d, k: __m256d) -> __m256d {
    _mm256_mul_pd(z, k)
}

/// The AVX2 form of `ComplexFft::radix8_stage`: current length `n` (at least 8), stride `s`,
/// with `s` a power of two and at least 2, so that two neighbouring butterflies `t` and `t + 1`
/// (`t` even) share `p` and the twiddle factors. Same tiles as the scalar stage, in pairs.
///
/// # Safety
/// The CPU must support AVX2.
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn radix8_stage(
    src: &[Complex],
    dst: &mut [Complex],
    twiddle: &[Complex],
    n: usize,
    s: usize,
    inverse: bool,
) {
    const TILE: usize = 8;
    let total = twiddle.len(); // the length of the plan
    assert!(s >= 2 && s.is_power_of_two() && n >= 8 && n * s == total);
    assert!(src.len() >= total && dst.len() >= total);
    let step = total / n;
    let span = total / 8; // distance of the input streams, and number of butterflies; even
    let shift = s.trailing_zeros();
    let len = TILE.min(span); // even, divides `span`
    let (sp, dp) = (src.as_ptr(), dst.as_mut_ptr());
    let r = _mm256_set1_pd(std::f64::consts::FRAC_1_SQRT_2);
    let zero = _mm256_setzero_pd();
    let mut a = [[zero; TILE / 2]; 4];
    let mut b = [[zero; TILE / 2]; 4];
    for t0 in (0..span).step_by(len) {
        for j in 0..4 {
            for u in 0..len / 2 {
                // Indices are below `8 * span = total`: `t0 + 2 * u + 1 < span`, `j + 4 < 8`.
                // SAFETY: in bounds by the assertion on the lengths, see above.
                let (x, y) = unsafe {
                    (
                        load(sp.add(j * span + t0 + 2 * u)),
                        load(sp.add((j + 4) * span + t0 + 2 * u)),
                    )
                };
                let d = _mm256_sub_pd(x, y);
                a[j][u] = _mm256_add_pd(x, y);
                // times exp(-+2 pi i j / 8): 1, (1 -+ i)/sqrt 2, -+i, (-1 -+ i)/sqrt 2
                b[j][u] = match j {
                    0 => d,
                    1 => scale(_mm256_add_pd(d, mul_j(d, inverse)), r),
                    2 => mul_j(d, inverse),
                    _ => scale(_mm256_sub_pd(mul_j(d, inverse), d), r),
                };
            }
        }
        for u in 0..len / 2 {
            let t = t0 + 2 * u;
            let (p, q) = (t >> shift, t & (s - 1)); // `q` is even and `q + 1 < s`
            let even = dft4(a[0][u], a[1][u], a[2][u], a[3][u], inverse);
            let odd = dft4(b[0][u], b[1][u], b[2][u], b[3][u], inverse);
            let out = q + s * 8 * p;
            for k in 0..4 {
                for (rr, v) in [(2 * k, even[k]), (2 * k + 1, odd[k])] {
                    let w = twiddle[rr * p * step];
                    let w = if inverse { w.conj() } else { w };
                    // SAFETY: `out + s * rr + 1 < s * 8 * (p + 1) <= s * n = total`.
                    unsafe {
                        store(
                            dp.add(out + s * rr),
                            mul(v, _mm256_set1_pd(w.re), _mm256_set1_pd(w.im)),
                        )
                    };
                }
            }
        }
    }
}
