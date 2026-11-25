# otfft-rs

A Rust port of [OTFFT](http://wwwa.pikara.ne.jp/okojisan/otfft-en/download.html) by OK Ojisan
(Takuya OKAHISA), based on the Stockham autosort FFT. MIT
license, with the notices of the original author and of the DEWETRON fork (see `LICENSE`).

The C++ fork and its measurements are the reference: <https://github.com/DEWETRON/otfft>
(`docs/benchmarks.md` has the speeds to compare with, `docs/spec/` the behaviour rules).

## Status

A scalar Stockham FFT for power-of-two lengths (classes `Dif4`, `Dif8`, `Dif16`, `Dit4`, `Dit8` and `Dit16`, step 3 in progress), with the names
and scales of the C++ library (`fwd` divides by N, `inv` does not). It is the correctness baseline,
not yet fast: 2.3x to 5.5x slower than RustFFT (`cargo bench`, `docs/benchmarks.md` B3).

```rust
use otfft_rs::{Complex, ComplexFft};

let fft = ComplexFft::new(1024).unwrap();
let mut x = vec![Complex::new(1.0, 0.0); 1024];
let mut scratch = vec![Complex::default(); fft.scratch_len()];
fft.process(&mut x, &mut scratch, otfft_rs::Direction::Forward, otfft_rs::Scale::Length).unwrap();
```

`cargo test` runs the checks against the direct DFT sum, round trips, all scales and the refused
inputs.

## Documentation

- [docs/spec/](docs/spec/README.md): what the crate does, rule by rule, with status tags and the tests
- [docs/decisions.md](docs/decisions.md): why it is built this way
- [docs/benchmarks.md](docs/benchmarks.md): measurements, with machine, compiler and command

## Plan

1. [x] Radix-2 Stockham, tests against the direct sum.
2. [x] Benchmark against RustFFT (criterion, as dev-dependencies only) and, with the optional
   feature `fftw` and the system `libfftw3`, against FFTW with the method of the C++ benchmark
   (`cargo bench --features fftw --bench vs_fftw`, B5, D11); `realfft` joins with step 5.
   Baseline: 4x to 14x slower than RustFFT ([docs/benchmarks.md](docs/benchmarks.md), B1).
3. [ ] Radix-4, 8 and 16 DIF/DIT kernels (the C++ classes DIF4 to DIT16), a table of the fastest
   class per length like `ffttune` produces. The best class differs between CPUs (C++
   benchmarks: Intel, Ryzen 9 9950X, Ryzen 9 3900X), so every class stays, even one that is slower
   on the development machine (D12), and step 8 chooses. Done: `Dif4`, `Dif8`, `Dif16`, `Dit4`, `Dit8` and `Dit16` (D8, D9, D12, D13, D14, D15, D16; on the
   Ryzen 3900X `Dif8` is fastest at 13 of 14 lengths, B6). The classes of step 3 are done; the table of the fastest class per length is step 8.
4. [ ] SIMD: AVX2 and SSE2 with run-time selection (`std::arch`, `is_x86_feature_detected!`),
   scalar fallback; later NEON.
5. [ ] Real FFT and DCT-II (C++: `RealFFT`, `RealDCT`).
6. [ ] Other lengths: mixed radix 3, 5, 7 and Bluestein for large prime factors (C++ roadmap R1.4,
   R1.5; keep the prime-length bug of issue #16 out: use `usize`/`u64` for index products).
7. [ ] Threads (rayon, optional feature), thresholds at the measured break-even (C++ D17).
8. [ ] Port of the `ffttune` tool (C++ `tools/ffttune`, R1.8, D12, D22, D23): a tool that times
   every kernel of step 3 per power-of-two length, checks each result against the direct sum or
   a reference kernel, prints a table, and writes a tuning file that the library can load to pick
   the fastest class per length. Needs the kernels of step 3 (the table of the fastest class
   there is first written by hand from the measurements); worth more after the SIMD kernels of
   step 4.

## Design notes

- A plan holds only the twiddle factors and is immutable, so it can be shared between threads. The
  scratch buffer comes from the caller (`process`); the `*_alloc` methods allocate one per call.
  The C++ library keeps the buffer in the object and says "one object per thread".
- `Complex` is `#[repr(C)]` with the layout of `OTFFT::complex_t`, so data can be exchanged with
  the C++ library without copying (for comparing results across the two).
- No dependencies in the library.
