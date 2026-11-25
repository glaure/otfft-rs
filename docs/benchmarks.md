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

B4: the last radix-4 stage without twiddles, and tiling it
-----------------------------------------------------------
*   **Machine, compiler, crates, command:** as B3 (`--save-baseline <name>`, one pass per variant
    in the same session, load average under 2). N = 2^11 and 2^17 were added to the benchmark.
*   **Variants:** `main` is the code of B3. `A` replaces the general radix-4 stage by a last-stage
    function without twiddle multiplications (the stage only runs with `n = 4`, where `p = 0` and
    all three twiddles are 1). `B` is A with tiles of 8, reading the four input streams before
    writing the four output streams (as the radix-8 stage of D9).
*   **Last stage** is the stage that ends the run: `r4` for log2(N) mod 3 = 2, `r2` for 1, none
    for 0. Only the `r4` rows are changed by A and B; the other rows show the noise.

| N         | last | main      | A         | B         | A / main | B / main | B / A |
|----------:|:----:|----------:|----------:|----------:|---------:|---------:|------:|
| 16        | r2   | 87.6 ns   | 68.9 ns   | 68.2 ns   | 0.79     | 0.78     | 0.99  |
| 64        | none | 346 ns    | 338 ns    | 330 ns    | 0.98     | 0.96     | 0.98  |
| 256       | r4   | 1.60 µs   | 1.51 µs   | 1.55 µs   | 0.95     | 0.97     | 1.02  |
| 1024      | r2   | 7.88 µs   | 7.79 µs   | 7.76 µs   | 0.99     | 0.99     | 1.00  |
| 2048      | r4   | 17.0 µs   | 15.9 µs   | 16.4 µs   | 0.93     | 0.96     | 1.03  |
| 4096      | none | 39.5 µs   | 39.0 µs   | 38.6 µs   | 0.99     | 0.98     | 0.99  |
| 16384     | r4   | 195 µs    | 181 µs    | 179 µs    | 0.93     | 0.92     | 0.99  |
| 65536     | r2   | 868 µs    | 852 µs    | 848 µs    | 0.98     | 0.98     | 1.00  |
| 131072    | r4   | 1.83 ms   | 1.76 ms   | 1.78 ms   | 0.96     | 0.98     | 1.01  |
| 262144    | none | 4.06 ms   | 4.01 ms   | 3.93 ms   | 0.99     | 0.97     | 0.98  |
| 1048576   | r4   | 22.8 ms   | 22.9 ms   | 23.6 ms   | 1.01     | 1.04     | 1.03  |

**Reading:** the rows that A and B do not change (r2, none) move by 1% to 4% between passes, and
N = 16 by 21%, so differences of a few percent are noise. A takes 4% to 7% less time than main at
N = 256, 2048, 16384 and 131072 (the unchanged rows move 1% to 4%), and nothing at N = 2^20. B is
no better than A (0.98 to 1.03 of A) and slower at 2^20. Tiling does not help a stage with four
streams, where it helped the radix-8 stage with eight (B3); that fits the cache-set explanation
for B3 but does not prove it, since B3 also gained at N = 16.

B5: against RustFFT and FFTW, with the method of the C++ benchmark
-------------------------------------------------------------------
*   **Machine, compiler:** as B1 (Ryzen 9 3900X, Linux 7.0.0-44-generic, rustc 1.100.0), idle
    (load average 0.3 after the run). **FFTW:** 3.3.10 from the Debian package, one thread.
    **RustFFT:** 6.4.1. **Code:** `otfft-rs` after D10 (radix-8 stages in tiles).
*   **Command:** `cargo bench --features fftw --bench vs_fftw -- --noplot --save-baseline fftw1`
    (needs `libfftw3-dev`; D11). criterion with 10 samples of 3 s per case, median shown.
*   **What is timed:** as in the C++ benchmark: one forward unnormalised complex transform with
    the plan built. `otfft-rs` works in place and the time includes copying the input into its
    work buffer; RustFFT and FFTW work out of place and do not copy. FFTW is planned with
    `FFTW_MEASURE` and `FFTW_ESTIMATE`, planning limited to 10 s. Input: random numbers in
    [-1, 1]. Every result was compared with RustFFT before it was timed (relative error below
    1e-9, no case was refused). MFLOPS are "equivalent", 5 N log2(N) divided by the time, as in
    the C++ benchmark.
*   **Differences from the C++ benchmark:** criterion statistics instead of 11 batches of 2 ms;
    data not aligned to 64 bytes (FFTW's buffers are, ours come from `Vec`); a different FFTW
    build (the C++ numbers use 3.3.10 built from the tarball with AVX2); no 12-thread run; powers
    of two only.

| N         | otfft-rs  | rustfft   | fftw-measure | fftw-estimate | otfft-rs / fftw-measure | MFLOPS otfft-rs | MFLOPS fftw-measure |
|----------:|----------:|----------:|-------------:|--------------:|-----------------------:|----------------:|--------------------:|
| 16        | 73.5 ns   | 13.9 ns   | 15.3 ns   | 15.4 ns   |  4.8x |   4354 |  20891 |
| 32        | 132.0 ns  | 27.8 ns   | 33.1 ns   | 32.8 ns   |  4.0x |   6060 |  24180 |
| 64        | 345.3 ns  | 64.0 ns   | 61.3 ns   | 61.4 ns   |  5.6x |   5561 |  31323 |
| 128       | 751.4 ns  | 138.3 ns  | 125.1 ns  | 124.7 ns  |  6.0x |   5962 |  35806 |
| 256       | 1.53 µs   | 298.0 ns  | 276.1 ns  | 274.6 ns  |  5.5x |   6690 |  37082 |
| 512       | 3.80 µs   | 751.8 ns  | 627.0 ns  | 610.0 ns  |  6.1x |   6062 |  36748 |
| 1024      | 7.92 µs   | 1.79 µs   | 1.49 µs   | 1.46 µs   |  5.3x |   6466 |  34393 |
| 2048      | 16.59 µs  | 4.06 µs   | 3.60 µs   | 3.66 µs   |  4.6x |   6789 |  31267 |
| 4096      | 40.55 µs  | 9.54 µs   | 9.96 µs   | 9.57 µs   |  4.1x |   6061 |  24678 |
| 8192      | 86.88 µs  | 21.43 µs  | 20.43 µs  | 21.10 µs  |  4.3x |   6129 |  26062 |
| 16384     | 185.01 µs | 47.89 µs  | 56.31 µs  | 58.85 µs  |  3.3x |   6199 |  20367 |
| 32768     | 430.29 µs | 105.85 µs | 124.20 µs | 122.09 µs |  3.5x |   5712 |  19787 |
| 65536     | 872.08 µs | 244.52 µs | 270.15 µs | 263.05 µs |  3.2x |   6012 |  19407 |
| 131072    | 1.80 ms   | 549.40 µs | 578.26 µs | 574.04 µs |  3.1x |   6190 |  19267 |
| 262144    | 4.22 ms   | 1.21 ms   | 1.32 ms   | 1.35 ms   |  3.2x |   5594 |  17868 |
| 524288    | 10.62 ms  | 3.45 ms   | 5.16 ms   | 5.15 ms   |  2.1x |   4690 |   9647 |
| 1048576   | 24.27 ms  | 8.76 ms   | 12.37 ms  | 12.34 ms  |  2.0x |   4321 |   8479 |
| 2097152   | 52.93 ms  | 19.88 ms  | 26.51 ms  | 26.68 ms  |  2.0x |   4160 |   8305 |
| 4194304   | 110.07 ms | 48.96 ms  | 83.80 ms  | 250.79 ms |  1.3x |   4192 |   5506 |

(Times: median of criterion's estimate; the ratio is time of `otfft-rs` over time of FFTW with
`FFTW_MEASURE`, so 3.3x means 3.3 times slower.)

**Reading:** `otfft-rs` is 1.3x to 6.1x slower than FFTW with `FFTW_MEASURE`: 4.0x to 6.1x for N
up to 2048, 3.1x to 4.3x from 4096 to 262144, and 2.0x to 2.1x from 2^19 to 2^21; at 2^22 the gap
is 1.3x, because FFTW's time rises there (84 ms against 49 ms for RustFFT). RustFFT takes up to
20% more time than FFTW for N from 512 to 2048, and equal time or up to 15% less from 4096 to
262144; from 2^19 it takes 1.3x to 1.7x less. `FFTW_ESTIMATE` is as fast as `MEASURE` everywhere
except N = 2^22 (251 ms against 84 ms). Equivalent MFLOPS of `otfft-rs` stay between 4100 and
6800 while FFTW's fall from about 37000 (N = 256) to 5500 (2^22). The flat MFLOPS up to 2^18
suggest that the scalar code is limited by instruction count and not by memory (not tested), which
is what the SIMD step (4) addresses.

The C++ `docs/benchmarks.md` reports OTFFT as fast as FFTW up to N = 64 and within 1.3x up to
N = 16384 on the same machine model. `otfft-rs` is not there yet: that is the gap that steps 3
(radix 16, the table of the fastest class) and 4 (SIMD) are to close. These numbers are from a
different session and FFTW build than the C++ ones, so they are not a head-to-head comparison.

B6: the three classes side by side (Dif4, Dif8, Dif16)
-------------------------------------------------------
*   **Machine, compiler, crates:** as B1 (Ryzen 9 3900X, rustc 1.100.0, rustfft 6.4.1). Load
    average about 1 before and after.
*   **Command:** `cargo bench --bench complex_fft -- --warm-up-time 2 --measurement-time 5
    --noplot --save-baseline kern3`: one pass, every class in the same session (D12).
*   **What:** as B1 (in place, no scaling, no copy), median of criterion's estimate. One pass
    only: differences under about 5% are within the noise of this machine (B3, B4).
*   **Classes:** `Dif4` is not tiled (B2). `Dif8` and `Dif16` work on tiles of 8 butterflies
    (D9, D13).

| N         | Dif4      | Dif8      | Dif16     | fastest | Dif16 / Dif8 | Dif4 / Dif8 | rustfft   |
|----------:|----------:|----------:|----------:|:-------:|-------------:|------------:|----------:|
| 16        | 147 ns    | 78 ns     | 96 ns     | Dif8  | 1.23 | 1.89 | 14 ns     |
| 32        | 340 ns    | 151 ns    | 154 ns    | Dif8  | 1.02 | 2.24 | 27 ns     |
| 64        | 694 ns    | 358 ns    | 309 ns    | Dif16 | 0.86 | 1.94 | 62 ns     |
| 256       | 3.09 µs   | 1.55 µs   | 1.81 µs   | Dif8  | 1.17 | 1.99 | 293 ns    |
| 1024      | 14.03 µs  | 8.18 µs   | 8.60 µs   | Dif8  | 1.05 | 1.72 | 1.81 µs   |
| 2048      | 29.73 µs  | 17.13 µs  | 20.15 µs  | Dif8  | 1.18 | 1.74 | 4.08 µs   |
| 4096      | 61.25 µs  | 41.70 µs  | 45.86 µs  | Dif8  | 1.10 | 1.47 | 9.53 µs   |
| 8192      | 131.52 µs | 87.01 µs  | 92.30 µs  | Dif8  | 1.06 | 1.51 | 21.94 µs  |
| 16384     | 279.74 µs | 188.98 µs | 193.55 µs | Dif8  | 1.02 | 1.48 | 50.10 µs  |
| 32768     | 589.53 µs | 432.74 µs | 449.65 µs | Dif8  | 1.04 | 1.36 | 106.30 µs |
| 65536     | 1.20 ms   | 891.76 µs | 997.65 µs | Dif8  | 1.12 | 1.35 | 237.10 µs |
| 131072    | 2.55 ms   | 1.84 ms   | 2.16 ms   | Dif8  | 1.17 | 1.39 | 525.60 µs |
| 262144    | 5.31 ms   | 4.25 ms   | 4.63 ms   | Dif8  | 1.09 | 1.25 | 1.14 ms   |
| 1048576   | 28.90 ms  | 24.27 ms  | 36.06 ms  | Dif8  | 1.49 | 1.19 | 9.81 ms   |

**Reading:** on this CPU `Dif8` is the fastest class at 13 of the 14 lengths. `Dif16` is faster
only at N = 64 (14% less time; one radix-16 stage and one radix-4 stage), within 5% of `Dif8` at
N = 32, 1024, 16384 and 32768, and 6% to 49% slower elsewhere; the largest loss is N = 2^20 (1.49x),
where every stage is radix 16. `Dif4` takes 1.2x to 2.2x the time of `Dif8`, with the gap
narrowing from 2.2x (N = 32) to 1.2x (2^20). Nothing here says what the order is on another CPU:
in the C++ benchmarks the order of `DIF4`, `DIF8` and `DIF16` changes between an Intel CPU and two
Ryzen CPUs (D12). That is the reason all three stay.


B7: the Dit4 class next to Dif4 and Dif8
----------------------------------------
*   **Machine, compiler, crates:** as B1 (Ryzen 9 3900X, rustc 1.100.0, rustfft 6.4.1). Load
    average 0.6 to 2 during the run (a single pass, so noise of about 5% applies, B3, B4).
*   **Command:** `cargo bench --bench complex_fft -- --noplot`: every class in the same session
    (D12), in place, no scaling, no copy, median of criterion's estimate.
*   **Classes:** `Dit4` is not tiled, like `Dif4` (D14).

| N         | Dif4      | Dif8      | Dit4      | Dit4 / Dif4 | Dit4 / Dif8 |
|----------:|----------:|----------:|----------:|------------:|------------:|
| 16        | 155 ns    | 74 ns     | 127 ns    | 0.82 | 1.73 |
| 32        | 334 ns    | 136 ns    | 320 ns    | 0.96 | 2.36 |
| 64        | 702 ns    | 330 ns    | 639 ns    | 0.91 | 1.93 |
| 256       | 3.12 µs   | 1.56 µs   | 2.92 µs   | 0.94 | 1.87 |
| 1024      | 13.97 µs  | 7.88 µs   | 13.20 µs  | 0.95 | 1.68 |
| 2048      | 30.57 µs  | 16.35 µs  | 28.36 µs  | 0.93 | 1.73 |
| 4096      | 61.56 µs  | 39.08 µs  | 57.88 µs  | 0.94 | 1.48 |
| 8192      | 135.5 µs  | 84.2 µs   | 130.6 µs  | 0.96 | 1.55 |
| 16384     | 284.7 µs  | 176.8 µs  | 270.3 µs  | 0.95 | 1.53 |
| 32768     | 603.8 µs  | 415.6 µs  | 560.2 µs  | 0.93 | 1.35 |
| 65536     | 1.189 ms  | 837 µs    | 1.142 ms  | 0.96 | 1.36 |
| 131072    | 2.547 ms  | 1.745 ms  | 2.507 ms  | 0.98 | 1.44 |
| 262144    | 5.194 ms  | 3.946 ms  | 5.094 ms  | 0.98 | 1.29 |
| 1048576   | 29.03 ms  | 22.88 ms  | 28.45 ms  | 0.98 | 1.24 |

**Reading:** `Dit4` takes 2% to 18% less time than `Dif4` at every length (the largest gain at
N = 16 and 64, within noise at 2^17 to 2^20), and 1.2x to 2.4x the time of `Dif8`, which stays the
fastest class on this CPU. The twiddle multiplication before the butterfly instead of after it
is the only difference, so the gain is small. Nothing here says what the order is on another CPU.


B8: the Dit8 class next to Dif8 and Dit4
----------------------------------------
*   **Machine, compiler, crates:** as B1 (Ryzen 9 3900X, rustc 1.100.0, rustfft 6.4.1). Load
    average 1 to 2 during the run (a single pass, so noise of about 5% applies, B3, B4).
*   **Command:** `cargo bench --bench complex_fft -- --noplot 'otfft-rs-(Dif8|Dit4|Dit8)'`, in
    place, no scaling, no copy, median of criterion's estimate. The `Dif8` and `Dit4` numbers are
    of this run, so they differ a little from B6 and B7.
*   **Classes:** `Dit8` is not tiled (D15).

| N         | Dif8      | Dit4      | Dit8      | Dit8 / Dif8 | Dit8 / Dit4 |
|----------:|----------:|----------:|----------:|------------:|------------:|
| 16        | 90 ns     | 133 ns    | 74 ns     | 0.81 | 0.55 |
| 32        | 136 ns    | 328 ns    | 147 ns    | 1.08 | 0.45 |
| 64        | 353 ns    | 646 ns    | 391 ns    | 1.11 | 0.61 |
| 256       | 1.57 µs   | 3.01 µs   | 1.78 µs   | 1.14 | 0.59 |
| 1024      | 8.14 µs   | 13.62 µs  | 9.12 µs   | 1.12 | 0.67 |
| 2048      | 16.60 µs  | 29.13 µs  | 20.00 µs  | 1.20 | 0.69 |
| 4096      | 40.99 µs  | 59.97 µs  | 57.46 µs  | 1.40 | 0.96 |
| 8192      | 88.19 µs  | 130.59 µs | 106.53 µs | 1.21 | 0.82 |
| 16384     | 181.02 µs | 276.28 µs | 225.73 µs | 1.25 | 0.82 |
| 32768     | 440.94 µs | 577.60 µs | 613.11 µs | 1.39 | 1.06 |
| 65536     | 884.48 µs | 1.16 ms   | 1.22 ms   | 1.38 | 1.05 |
| 131072    | 1.82 ms   | 2.58 ms   | 2.52 ms   | 1.38 | 0.98 |
| 262144    | 4.12 ms   | 5.30 ms   | 6.26 ms   | 1.52 | 1.18 |
| 1048576   | 23.45 ms  | 30.33 ms  | 28.09 ms  | 1.20 | 0.93 |

**Reading:** `Dit8` takes 0.55 to 0.96 of the time of `Dit4` up to N = 2^14 (2^12 is 0.96, within
noise) and 0.93 to 1.18 above it. Against `Dif8` it wins only at N = 16 (0.81) and takes 1.08x to
1.52x the time elsewhere, so `Dif8` stays the fastest class on this CPU. Nothing here says what
the order is on another CPU.
