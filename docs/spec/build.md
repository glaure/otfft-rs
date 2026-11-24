Build and platforms
===================

Rules
-----

*   **BLD-1** `[done]` `cargo test` builds the library and runs the unit tests, the integration
    tests in `tests/` and the doctest in `src/lib.rs`.
*   **BLD-2** `[built]` The library has no dependencies. `[dev-dependencies]` stay empty until the
    benchmark step adds `criterion`, `rustfft` and `realfft` there (README, step 2).
*   **BLD-3** `[built]` Edition 2024 (Rust 1.85 or newer). Developed with Rust 1.100 on Windows 11.
*   **BLD-4** `[done]` `publish = false` in `Cargo.toml`: `cargo publish` refuses, the crate is private.
*   **BLD-5** `[built]` `cargo clippy` reports nothing on the current code.
*   **BLD-6** `[gap]` No CI yet. Right: build and test on Linux and Windows.

Decided
-------

*   D4: private crate.

Questions
---------

*   Where is the repository hosted, and which CI does it use?
