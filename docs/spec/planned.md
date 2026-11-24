Not ported yet
==============

Areas of the C++ library (spec in the DEWETRON repository, `docs/spec/`) and what is planned for
them here. Order as in the README plan.

Rules
-----

*   **PLN-1** `[planned]` Radix-4, 8 and 16 kernels, DIF and DIT, and a table of the fastest class per
    length (C++: DIF4 to DIT16, `tools/ffttune`, D12, D14). README step 3.
*   **PLN-2** `[planned]` SIMD with run-time selection: AVX2 and SSE2 with `std::arch` and
    `is_x86_feature_detected!`, a scalar fallback, NEON later (C++ D3). README step 4.
*   **PLN-3** `[planned]` Real FFT and DCT-II (C++: `RealFFT` with N output values, input of `inv`
    destroyed; `RealDCT`). README step 5.
*   **PLN-4** `[planned]` Any length: mixed radix with 3, 5 and 7, and Bluestein for a prime factor
    above 100 (C++ CFFT-5, CFFT-7, D13). Index products use `u64`, so the overflow of C++ issue #16
    cannot happen. README step 6.
*   **PLN-5** `[planned]` Threads through an optional `rayon` feature, with thresholds at the
    measured break-even (C++ CFFT-11, D17). README step 7.
*   **PLN-6** `[planned]` Benchmarks against RustFFT and `realfft` with `criterion`, and a table
    that compares with the C++ numbers in `docs/benchmarks.md` of the fork. README step 2.

Questions
---------

*   Does the port aim to match the C++ results bit for bit (same operation order) or only within
    rounding? Today only within rounding is tested.
