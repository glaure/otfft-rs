Build and platforms
===================

Rules
-----

*   **BLD-1** `[done]` `cargo test` builds the library and runs the unit tests, the integration
    tests in `tests/` and the doctest in `src/lib.rs`.
*   **BLD-2** `[built]` The library has no dependencies. `[dev-dependencies]` hold `criterion`,
    `rustfft` and `realfft` for the benchmarks (D7); `realfft` is not used before the real FFT.
*   **BLD-3** `[built]` Edition 2024 (Rust 1.85 or newer). Developed with Rust 1.100 on Windows 11.
*   **BLD-4** `[done]` `publish = false` in `Cargo.toml`: `cargo publish` refuses, the crate is private.
*   **BLD-5** `[built]` `cargo clippy` reports nothing on the current code.
*   **BLD-6** `[built]` CI (GitHub Actions) runs `cargo test` and `cargo clippy --all-targets` on
    Linux, Windows and macOS. It does not run benchmarks (D6).
*   **BLD-7** `[built]` `cargo bench --bench complex_fft` times the complex forward transform
    against RustFFT. The results are in `docs/benchmarks.md` (D7).
*   **BLD-8** `[done]` The cargo feature `fftw` (off by default) links the system `libfftw3` for
    `cargo bench --features fftw --bench vs_fftw`, which compares `otfft-rs`, RustFFT and FFTW
    (`FFTW_MEASURE` and `FFTW_ESTIMATE`) on power-of-two lengths 2^4 to 2^22 with the method of the
    C++ benchmark, and for `cargo test --features fftw`. The library does not use FFTW; without
    the feature nothing links it and CI does not need it (D11). (`the_binding_computes_the_direct_sum`,
    `otfft_rs_equals_fftw` in `tests/fftw.rs`)

Decided
-------

*   D4: private crate.
*   D6: CI on three operating systems.
*   D7: benchmark against RustFFT with criterion.
*   D11: benchmark against FFTW through an optional feature.

Questions
---------

*   None open.
