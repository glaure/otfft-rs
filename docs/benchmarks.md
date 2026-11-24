Benchmarks
==========

Rule: a speed claim in a doc or commit message needs a line here (AGENTS.md rule 9). Compare with
the C++ library only on the same machine.

B1: complex forward, otfft-rs against RustFFT (step 0 baseline)
----------------------------------------------------------------
*   **Machine:** AMD Ryzen 9 3900X (12 cores, AVX2), Linux 7.0.0-44-generic.
*   **Compiler:** rustc 1.100.0 (0ef1bc386 2026-11-09), `bench` profile (release), no
    `target-cpu` flag.
*   **Crates:** rustfft 6.4.1, criterion 0.8.2.
*   **Command:** `cargo bench --bench complex_fft -- --warm-up-time 1 --measurement-time 2`
*   **What:** `f64` complex, in place, forward, no scaling, caller-owned scratch, power-of-two
    lengths. Median of criterion's estimate. Short runs (2 s): treat differences under ~10% as
    noise.

| N         | otfft-rs  | rustfft   | otfft-rs / rustfft |
|----------:|----------:|----------:|-------------------:|
| 16        | 189 ns    | 14.6 ns   | 12.9x slower       |
| 64        | 858 ns    | 63.1 ns   | 13.6x              |
| 256       | 3.80 µs   | 300 ns    | 12.7x              |
| 1024      | 17.3 µs   | 2.20 µs   | 7.9x               |
| 4096      | 99.6 µs   | 13.8 µs   | 7.2x               |
| 16384     | 392 µs    | 71.7 µs   | 5.5x               |
| 65536     | 1.72 ms   | 338 µs    | 5.1x               |
| 262144    | 6.34 ms   | 1.38 ms   | 4.6x               |
| 1048576   | 60.9 ms   | 15.2 ms   | 4.0x               |

**Reading:** the scalar radix-2 baseline is 4x to 14x slower than RustFFT. The gap is largest for
small N, where radix-2 does the most passes and RustFFT uses unrolled butterflies. No C++ numbers
yet: they need the C++ library built and run on this machine.

B2: radix-4 stages (step 3, first kernel)
-----------------------------------------
*   **Machine, compiler, crates:** as B1 (Ryzen 9 3900X, rustc 1.100.0, rustfft 6.4.1).
*   **Command:** `cargo bench --bench complex_fft -- --warm-up-time 2 --measurement-time 6 --noplot`
*   **What:** as B1, with radix-4 Stockham stages and one radix-2 stage when log2(N) is odd (D8).
    The runs of B1 and B2 are different sessions on a shared machine: RustFFT's own time moved by
    up to 25% between two B2 trial runs. The ratio to RustFFT in the same run is the number to
    read; N = 64 and 256 were disturbed in this run (both libraries 17% to 29% slower than in the
    previous trial), so treat those two rows as +-15%.

| N         | otfft-rs  | rustfft   | otfft-rs / rustfft | B1 ratio |
|----------:|----------:|----------:|-------------------:|---------:|
| 16        | 168 ns    | 15.0 ns   | 11.2x slower       | 12.9x    |
| 64        | 876 ns    | 73.9 ns   | 11.9x              | 13.6x    |
| 256       | 3.82 µs   | 311 ns    | 12.3x              | 12.7x    |
| 1024      | 14.8 µs   | 1.84 µs   | 8.1x               | 7.9x     |
| 4096      | 63.9 µs   | 9.74 µs   | 6.6x               | 7.2x     |
| 16384     | 289 µs    | 49.3 µs   | 5.9x               | 5.5x     |
| 65536     | 1.21 ms   | 248 µs    | 4.9x               | 5.1x     |
| 262144    | 5.42 ms   | 1.21 ms   | 4.5x               | 4.6x     |
| 1048576   | 29.6 ms   | 9.94 ms   | 3.0x               | 4.0x     |

**Reading:** the absolute time of otfft-rs falls by about 15% to 50% (1024: 17.3 to 14.8 µs;
2^20: 60.9 to 29.6 ms), but the ratio to RustFFT hardly moves up to N = 65536, which is inside
the noise. Only N = 2^20 clearly gains on RustFFT (4.0x to 3.0x). Radix 4 alone does not close
the gap; RustFFT's lead at small N comes from unrolled butterflies, which radix 8 and 16 and SIMD
(steps 3 and 4) address.
