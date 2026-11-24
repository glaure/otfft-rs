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

B3: radix-8 stages (step 3, second kernel)
------------------------------------------
*   **Machine, compiler, crates:** as B1 (Ryzen 9 3900X, rustc 1.100.0, rustfft 6.4.1).
*   **Command:** `cargo bench --bench complex_fft -- --warm-up-time 2 --measurement-time 5 --noplot
    --save-baseline <name>`, once per variant, one after the other in the same session on an idle
    machine (load average under 2, RustFFT within 5% between the passes, except N = 16 at 17%).
*   **Variants:** `radix4` is the code of B2 (`main` at 19f5fbc). `radix8 v1` is a direct
    radix-8 stage (one butterfly per step, arrays of 8 values). `radix8 tiled` is the stage of D9
    (tiles of 8 elements, a pair of input streams at a time).

| N         | radix4    | radix8 v1 | radix8 tiled | tiled / radix4 | tiled / rustfft | radix4 / rustfft |
|----------:|----------:|----------:|-------------:|---------------:|----------------:|-----------------:|
| 16        | 165 ns    | 166 ns    | 84.0 ns      | 0.51           | 5.4x            | 12.4x            |
| 64        | 729 ns    | 779 ns    | 337 ns       | 0.46           | 5.3x            | 11.7x            |
| 256       | 3.24 µs   | 3.26 µs   | 1.56 µs      | 0.48           | 5.5x            | 10.7x            |
| 1024      | 14.6 µs   | 14.7 µs   | 7.83 µs      | 0.54           | 4.3x            | 8.1x             |
| 4096      | 62.4 µs   | 70.4 µs   | 39.2 µs      | 0.63           | 4.0x            | 6.3x             |
| 16384     | 288 µs    | 313 µs    | 187 µs       | 0.65           | 3.8x            | 5.7x             |
| 65536     | 1.23 ms   | 1.44 ms   | 847 µs       | 0.69           | 3.4x            | 5.1x             |
| 262144    | 5.41 ms   | 6.55 ms   | 3.95 ms      | 0.73           | 3.3x            | 4.5x             |
| 1048576   | 29.4 ms   | 31.9 ms   | 22.7 ms      | 0.77           | 2.3x            | 2.9x             |

The ratios to RustFFT use the RustFFT time of the same pass.

**Reading:** the direct radix-8 stage (v1) was no faster than radix 4 and 9% to 21% slower from
N = 4096. The tiled stage takes 23% to 54% less time than radix 4 at every length, and is 2.3x to
5.5x slower than RustFFT. The gain also appears at N = 16, where cache sets cannot matter, so the
tiling does more than avoid cache-set conflicts; the cause of the slow v1 was not isolated
(candidates: the streams `N/8` apart share cache sets, the arrays and closures of v1, and bounds
checks). The same tiling may help the radix-4 stage; not tried.
