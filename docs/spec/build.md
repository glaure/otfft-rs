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

Decided
-------

*   D4: private crate.
*   D6: CI on three operating systems.
*   D7: benchmark against RustFFT with criterion.

Questions
---------

*   None open.
