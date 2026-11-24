//! Stockham FFT after OTFFT by OK Ojisan (Takuya OKAHISA).
//!
//! Scalar Stockham autosort FFT for power-of-two lengths: radix-8 stages, then one radix-4 or
//! radix-2 stage when `log2(N)` is not a multiple of 3 (step 3). It is the correctness baseline
//! that the radix-16 and SIMD versions are measured against.
//!
//! The names and scales follow the C++ library (`fwd` divides by N, `inv` does not):
//!
//! | forward | inverse | scale |
//! | --- | --- | --- |
//! | `fwd`, `fwdn` | `invn` | 1/N |
//! | `fwd0` | `inv`, `inv0` | none |
//! | `fwdu` | `invu` | 1/sqrt(N) |
//!
//! ```
//! use otfft_rs::{Complex, ComplexFft};
//!
//! let fft = ComplexFft::new(8).unwrap();
//! let mut x = vec![Complex::new(1.0, 0.0); 8];
//! fft.fwd_alloc(&mut x); // a constant signal: all of the energy is in bin 0
//! assert!((x[0].re - 1.0).abs() < 1e-12);
//! assert!(x[1].abs() < 1e-12);
//! fft.inv_alloc(&mut x); // inv undoes fwd
//! assert!((x[3].re - 1.0).abs() < 1e-12);
//! ```

use std::f64::consts::PI;
use std::fmt;
use std::ops::{Add, Mul, Sub};

/// A complex number of two `f64`, with the same layout as `OTFFT::complex_t` of the C++ library.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(C)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub const fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }

    pub fn conj(self) -> Self {
        Complex::new(self.re, -self.im)
    }

    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// `exp(i * theta)`.
    pub fn expj(theta: f64) -> Self {
        let (s, c) = theta.sin_cos();
        Complex::new(c, s)
    }

    fn scale(self, k: f64) -> Self {
        Complex::new(self.re * k, self.im * k)
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, o: Complex) -> Complex {
        Complex::new(self.re + o.re, self.im + o.im)
    }
}

impl Sub for Complex {
    type Output = Complex;
    fn sub(self, o: Complex) -> Complex {
        Complex::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for Complex {
    type Output = Complex;
    fn mul(self, o: Complex) -> Complex {
        Complex::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}

/// Why a plan or a transform call was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    /// The length is 0 or not a power of two (other lengths come with the mixed radix and Bluestein steps).
    UnsupportedLength(usize),
    /// A buffer does not have the length of the plan.
    WrongBufferLength { expected: usize, got: usize },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::UnsupportedLength(n) => write!(f, "length {n} is not a power of two"),
            Error::WrongBufferLength { expected, got } => {
                write!(f, "buffer has {got} elements, the plan is for {expected}")
            }
        }
    }
}

impl std::error::Error for Error {}

/// Direction of the transform: the sign of the exponent.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Direction {
    Forward,
    Inverse,
}

/// Scale applied to the result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scale {
    None,
    /// 1/N
    Length,
    /// 1/sqrt(N)
    Unitary,
}

/// Plan of a complex in-place FFT of a power-of-two length. The plan holds only the twiddle
/// factors, so it can be shared between threads; every call brings its own scratch buffer.
#[derive(Clone, Debug)]
pub struct ComplexFft {
    n: usize,
    /// `exp(-2 pi i k / n)` for `k` in `0..n` (the radix-4 stage reads up to index `3n/4`)
    twiddle: Vec<Complex>,
}

impl ComplexFft {
    /// Plan for length `n` (a power of two, at least 1).
    pub fn new(n: usize) -> Result<Self, Error> {
        if n == 0 || !n.is_power_of_two() {
            return Err(Error::UnsupportedLength(n));
        }
        let twiddle = (0..n).map(|k| Complex::expj(-2.0 * PI * k as f64 / n as f64)).collect();
        Ok(ComplexFft { n, twiddle })
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        false
    }

    /// Length of the scratch buffer that `process` needs.
    pub fn scratch_len(&self) -> usize {
        self.n
    }

    /// In-place transform of `x` with `scratch` as work space (both of length `len()`).
    pub fn process(
        &self,
        x: &mut [Complex],
        scratch: &mut [Complex],
        dir: Direction,
        scale: Scale,
    ) -> Result<(), Error> {
        for buf in [&*x, &*scratch] {
            if buf.len() != self.n {
                return Err(Error::WrongBufferLength { expected: self.n, got: buf.len() });
            }
        }
        self.stockham(x, scratch, dir);
        let k = match scale {
            Scale::None => return Ok(()),
            Scale::Length => 1.0 / self.n as f64,
            Scale::Unitary => 1.0 / (self.n as f64).sqrt(),
        };
        for v in x.iter_mut() {
            *v = v.scale(k);
        }
        Ok(())
    }

    /// Stockham autosort, decimation in frequency. Each radix-8 stage splits the transform of
    /// length `n` with stride `s` into eight of length `n/8` with stride `8s`; one radix-4 or
    /// radix-2 stage ends the run when `log2(N)` is not a multiple of 3. The data moves between
    /// `x` and `y`, and no bit reversal is needed.
    fn stockham(&self, x: &mut [Complex], y: &mut [Complex], dir: Direction) {
        let total = self.n;
        let inverse = dir == Direction::Inverse;
        let mut n = total;
        let mut s = 1;
        let mut in_x = true; // where the current data is
        while n > 1 {
            let (src, dst): (&[Complex], &mut [Complex]) =
                if in_x { (&*x, &mut *y) } else { (&*y, &mut *x) };
            if n >= 8 {
                self.radix8_stage(src, dst, n, s, inverse);
                n /= 8;
                s *= 8;
            } else if n == 4 {
                self.radix4_last(src, dst, s, inverse);
                n = 1;
            } else {
                // n == 2: the only twiddle factor is 1
                for q in 0..s {
                    let a = src[q];
                    let b = src[q + s];
                    dst[q] = a + b;
                    dst[q + s] = a - b;
                }
                n = 1;
            }
            in_x = !in_x;
        }
        if !in_x {
            x.copy_from_slice(y);
        }
    }

    /// One radix-8 DIF stage of the Stockham FFT: `src` to `dst`, current length `n`, stride `s`.
    /// The eight inputs are folded into two radix-4 butterflies: the sums `x[j] + x[j+4]` give the
    /// even outputs, the differences times `exp(-+2 pi i j / 8)` give the odd ones.
    ///
    /// The eight input streams are `N/8` elements apart, a power of two, so they fall in the same
    /// cache sets. The stage therefore works on tiles of `TILE` elements: it reads a pair of
    /// streams at a time, and finishes a cache line before it goes to the next stream.
    fn radix8_stage(&self, src: &[Complex], dst: &mut [Complex], n: usize, s: usize, inverse: bool) {
        const TILE: usize = 8;
        let step = self.n / n;
        let span = self.n / 8; // distance of the input streams, and number of butterflies
        let shift = s.trailing_zeros(); // `s` is a power of two: p = t >> shift, q = t & (s - 1)
        let len = TILE.min(span);
        // multiplication by -i (forward) or +i (inverse)
        let mj = |z: Complex| {
            if inverse { Complex::new(-z.im, z.re) } else { Complex::new(z.im, -z.re) }
        };
        let r = std::f64::consts::FRAC_1_SQRT_2;
        let dft4 = |u: [Complex; 4]| {
            let (e, f) = (u[0] + u[2], u[0] - u[2]);
            let (g, h) = (u[1] + u[3], mj(u[1] - u[3]));
            [e + g, f + h, e - g, f - h]
        };
        let mut a = [[Complex::default(); TILE]; 4];
        let mut b = [[Complex::default(); TILE]; 4];
        for t0 in (0..span).step_by(len) {
            for j in 0..4 {
                let lo = &src[j * span + t0..j * span + t0 + len];
                let hi = &src[(j + 4) * span + t0..(j + 4) * span + t0 + len];
                for u in 0..len {
                    let (x, y) = (lo[u], hi[u]);
                    let d = x - y;
                    a[j][u] = x + y;
                    // times exp(-+2 pi i j / 8): 1, (1 -+ i)/sqrt 2, -+i, (-1 -+ i)/sqrt 2
                    b[j][u] = match j {
                        0 => d,
                        1 => (d + mj(d)).scale(r),
                        2 => mj(d),
                        _ => (mj(d) - d).scale(r),
                    };
                }
            }
            for u in 0..len {
                let t = t0 + u;
                let (p, q) = (t >> shift, t & (s - 1));
                let (even, odd) = (dft4([a[0][u], a[1][u], a[2][u], a[3][u]]),
                                   dft4([b[0][u], b[1][u], b[2][u], b[3][u]]));
                let out = q + s * 8 * p;
                for k in 0..4 {
                    for (r, v) in [(2 * k, even[k]), (2 * k + 1, odd[k])] {
                        let w = self.twiddle[r * p * step];
                        let w = if inverse { w.conj() } else { w };
                        dst[out + s * r] = v * w;
                    }
                }
            }
        }
    }

    /// The last stage when `log2(N)` leaves a remainder of 2: radix 4 with `n = 4`, so `p = 0` and
    /// all three twiddle factors are 1. `s = N/4` is the stride.
    fn radix4_last(&self, src: &[Complex], dst: &mut [Complex], s: usize, inverse: bool) {
        for q in 0..s {
            let a = src[q];
            let b = src[q + s];
            let c = src[q + 2 * s];
            let d = src[q + 3 * s];
            let apc = a + c;
            let amc = a - c;
            let bpd = b + d;
            let bmd = b - d;
            // multiplication of `bmd` by -i (forward) or +i (inverse)
            let jbmd = if inverse {
                Complex::new(-bmd.im, bmd.re)
            } else {
                Complex::new(bmd.im, -bmd.re)
            };
            dst[q] = apc + bpd;
            dst[q + s] = amc + jbmd;
            dst[q + 2 * s] = apc - bpd;
            dst[q + 3 * s] = amc - jbmd;
        }
    }

    fn run_alloc(&self, x: &mut [Complex], dir: Direction, scale: Scale) {
        let mut scratch = vec![Complex::default(); self.n];
        self.process(x, &mut scratch, dir, scale).expect("buffer length");
    }

    // The names of the C++ library. The `_alloc` versions allocate the scratch buffer on each
    // call; use `process` to reuse one.

    /// Forward, scaled by 1/N (like `OTFFT::ComplexFFT::fwd`).
    pub fn fwd_alloc(&self, x: &mut [Complex]) {
        self.run_alloc(x, Direction::Forward, Scale::Length);
    }
    /// Forward, no scale.
    pub fn fwd0_alloc(&self, x: &mut [Complex]) {
        self.run_alloc(x, Direction::Forward, Scale::None);
    }
    /// Forward, scaled by 1/sqrt(N).
    pub fn fwdu_alloc(&self, x: &mut [Complex]) {
        self.run_alloc(x, Direction::Forward, Scale::Unitary);
    }
    /// Inverse, no scale (so `inv(fwd(x)) == x`).
    pub fn inv_alloc(&self, x: &mut [Complex]) {
        self.run_alloc(x, Direction::Inverse, Scale::None);
    }
    /// Inverse, scaled by 1/sqrt(N).
    pub fn invu_alloc(&self, x: &mut [Complex]) {
        self.run_alloc(x, Direction::Inverse, Scale::Unitary);
    }
    /// Inverse, scaled by 1/N.
    pub fn invn_alloc(&self, x: &mut [Complex]) {
        self.run_alloc(x, Direction::Inverse, Scale::Length);
    }
}
