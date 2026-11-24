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
