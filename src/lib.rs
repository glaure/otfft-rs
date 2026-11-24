//! Stockham FFT after OTFFT by OK Ojisan (Takuya OKAHISA).
//!
//! Scalar Stockham autosort FFT for power-of-two lengths: radix-4 stages (step 3, first kernel)
//! and one radix-2 stage when `log2(N)` is odd. It is the correctness baseline that the radix-8/16
//! and SIMD versions are measured against.
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

    /// Stockham autosort, decimation in frequency. Each radix-4 stage splits the transform of
    /// length `n` with stride `s` into four of length `n/4` with stride `4s`; when `log2(N)` is
    /// odd, one radix-2 stage ends the run. The data moves between `x` and `y`, and no bit
    /// reversal is needed.
    fn stockham(&self, x: &mut [Complex], y: &mut [Complex], dir: Direction) {
        let total = self.n;
        let inverse = dir == Direction::Inverse;
        let mut n = total;
        let mut s = 1;
        let mut in_x = true; // where the current data is
        while n > 1 {
            let (src, dst): (&[Complex], &mut [Complex]) =
                if in_x { (&*x, &mut *y) } else { (&*y, &mut *x) };
            if n >= 4 {
                self.radix4_stage(src, dst, n, s, inverse);
                n /= 4;
                s *= 4;
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

    /// One radix-4 DIF stage of the Stockham FFT: `src` to `dst`, current length `n`, stride `s`.
    fn radix4_stage(&self, src: &[Complex], dst: &mut [Complex], n: usize, s: usize, inverse: bool) {
        let m = n / 4;
        let step = self.n / n; // twiddle index of `exp(-2 pi i / n)`
        for p in 0..m {
            let tw = |k: usize| {
                let w = self.twiddle[k * step];
                if inverse { w.conj() } else { w }
            };
            let (w1, w2, w3) = (tw(p), tw(2 * p), tw(3 * p));
            for q in 0..s {
                let a = src[q + s * p];
                let b = src[q + s * (p + m)];
                let c = src[q + s * (p + 2 * m)];
                let d = src[q + s * (p + 3 * m)];
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
                dst[q + s * (4 * p)] = apc + bpd;
                dst[q + s * (4 * p + 1)] = (amc + jbmd) * w1;
                dst[q + s * (4 * p + 2)] = (apc - bpd) * w2;
                dst[q + s * (4 * p + 3)] = (amc - jbmd) * w3;
            }
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
