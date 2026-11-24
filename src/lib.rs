//! Stockham FFT after OTFFT by OK Ojisan (Takuya OKAHISA).
//!
//! Step 0: scalar radix-2 Stockham autosort FFT for power-of-two lengths. It is the
//! correctness baseline that the faster radix-4/8/16 and SIMD versions are measured against.
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
    /// `exp(-2 pi i k / n)` for `k` in `0..n/2`
    twiddle: Vec<Complex>,
}

impl ComplexFft {
    /// Plan for length `n` (a power of two, at least 1).
    pub fn new(n: usize) -> Result<Self, Error> {
        if n == 0 || !n.is_power_of_two() {
            return Err(Error::UnsupportedLength(n));
        }
        let twiddle = (0..n / 2).map(|k| Complex::expj(-2.0 * PI * k as f64 / n as f64)).collect();
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

    /// Stockham autosort, radix 2 decimation in frequency. At each stage the transform of
    /// length `n` with stride `s` is split into two of length `n/2` with stride `2s`, the data
    /// moves between `x` and `y`, and no bit reversal is needed.
    fn stockham(&self, x: &mut [Complex], y: &mut [Complex], dir: Direction) {
        let total = self.n;
        let mut n = total;
        let mut s = 1;
        let mut in_x = true; // where the current data is
        while n > 1 {
            let m = n / 2;
            let (src, dst): (&[Complex], &mut [Complex]) =
                if in_x { (&*x, &mut *y) } else { (&*y, &mut *x) };
            for p in 0..m {
                let w = self.twiddle[p * s];
                let w = if dir == Direction::Inverse { w.conj() } else { w };
                for q in 0..s {
                    let a = src[q + s * p];
                    let b = src[q + s * (p + m)];
                    dst[q + s * (2 * p)] = a + b;
                    dst[q + s * (2 * p + 1)] = (a - b) * w;
                }
            }
            n = m;
            s *= 2;
            in_x = !in_x;
        }
        if !in_x {
            x.copy_from_slice(y);
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
