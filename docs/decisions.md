Decisions
=========

Why `otfft-rs` is built the way it is. One entry per choice, newest at the end. A replaced entry
is marked superseded, not deleted.

D1: a port of the DEWETRON fork, read from the source
-----------------------------------------------------
*   **What:** `otfft-rs` ports the algorithms of OTFFT (OK Ojisan) as they are in the DEWETRON fork
    (<https://github.com/DEWETRON/otfft>), and uses its spec and benchmarks as the reference.
    Behaviour that differs is a decision here.
*   **Why:** the fork has tests, a spec and measured numbers for every area; a port can be checked
    against them.
*   **Rejected:** a binding to the C++ library (does not give a Rust-only build); a port from the
    original upstream only (no tests to compare with).

D2: C++ names and scales
------------------------
*   **What:** `fwd` divides by N and `inv` does not; `fwd0`, `fwdu`, `invu` and `invn` as in the
    C++ library. `process` takes `Direction` and `Scale` as arguments and the named methods are
    shortcuts.
*   **Why:** results can be compared with the C++ library and with its documentation without
    translating scales.
*   **Rejected:** RustFFT's convention (no scaling at all; the caller scales).

D3: immutable plan, scratch from the caller
-------------------------------------------
*   **What:** `ComplexFft` holds only the twiddle factors; `process(x, scratch, ...)` takes the
    scratch buffer. The plan is `Send + Sync`. `*_alloc` methods allocate a buffer per call.
*   **Why:** the C++ object owns a work buffer, so "one object per thread" is a rule the caller must
    remember (C++ FAC-6). In Rust the type system can allow sharing, and the caller can reuse one
    buffer for many calls.
*   **Rejected:** an internal buffer in a `RefCell` or a `Mutex` (not `Sync`, or a lock per call).
*   **Revisit when:** the threaded transforms (README step 7) need per-thread buffers inside a call.

D4: a private crate
-------------------
*   **What:** `publish = false`; no registry release. License MIT, with the notices of the original
    author and the fork kept.
*   **Why:** the project is private. The license of OTFFT (MIT) requires the notices to stay.

D5: scalar radix 2 first (stages replaced by D8)
------------------------------------------------
*   **What:** step 0 is a plain radix-2 Stockham with a twiddle table, power-of-two lengths only.
*   **Why:** a correct and tested baseline to measure the radix-4/8/16 and SIMD steps against, and a
    harness (direct-sum tests) that every later kernel must pass.
*   **Limit:** slow; not a replacement for the C++ library yet (CFFT-11).

D6: CI on three operating systems
---------------------------------
*   **What:** GitHub Actions runs `cargo test` and `cargo clippy --all-targets` on Linux, Windows and
    macOS with the stable toolchain; warnings are errors (`RUSTFLAGS=-D warnings`).
*   **Why:** the checks of AGENTS.md rule 10 run on every push and pull request, on every target
    the SIMD steps (README step 4) will have to support.
*   **Limit:** no benchmarks in CI; shared runners are too noisy for rule 9 measurements.

D7: benchmark against RustFFT with criterion
--------------------------------------------
*   **What:** `benches/complex_fft.rs` (criterion) times the complex forward transform of `otfft-rs`
    against RustFFT, in place with a caller-owned scratch buffer and no scaling, for N = 2^4 to
    2^20. `criterion`, `rustfft` and `realfft` are `[dev-dependencies]` only. Results go in
    `docs/benchmarks.md` with machine, compiler and command.
*   **Why:** RustFFT is the Rust reference for speed, and the numbers show which of README steps 3
    and 4 pays off first. Timing without allocation and scaling compares the transforms only.
*   **Rejected:** comparing only with the C++ library (needs it built on the same machine, and gives
    no Rust baseline); including scaling in the timing (RustFFT does not scale).
*   **Limit:** `realfft` is unused until the real FFT (README step 5). Short runs are noisy.

D8: radix-4 stages, with a radix-2 stage for odd log2(N) (stage order replaced by D9)
---------------------------------------------------------------------------------------
*   **What:** the transform runs radix-4 Stockham DIF stages (C++ class DIF4, written in scalar
    form) and ends with one radix-2 stage when log2(N) is odd. The twiddle table has N entries
    (it was N/2) because the radix-4 stage reads `exp(-2 pi i * 3p / n)`, up to index 3N/4. The
    radix-2 stage is no longer a general stage: it only runs with n = 2, where the twiddle is 1.
*   **Why:** README step 3, first kernel: half the passes over the data. Measured in B2: 15% to
    50% less time in absolute terms; not yet a clear gain over RustFFT below N = 2^20.
*   **Rejected:** a table of N/2 entries with a sign flip for indices above N/2 (a branch per
    twiddle for 4 KiB to 8 MiB saved); a runtime choice between radix 2 and 4 (no measured case
    where radix 2 wins; the per-length table of step 3 comes with radix 8 and 16).
*   **Limit:** the C++ `DIF4` is an AVX implementation; this is its scalar form, not a port of the
    SIMD code (step 4). Radix 8, 16 and DIT variants are separate commits.

D9: radix-8 stages in tiles, radix 4 or 2 to finish
----------------------------------------------------
*   **What:** the transform runs radix-8 Stockham DIF stages (C++ class DIF8, scalar form) while
    the current length is at least 8, then one radix-4 stage (length 4) or radix-2 stage (length
    2). The stage works on tiles of 8 butterflies: for each of the four pairs of streams `j` and
    `j + 4` it reads the tile, forms the sum and the twiddled difference, and only then does the
    two radix-4 butterflies and the stores. The radix-4 stage of D8 remains for the last stage (D10).
*   **Why:** README step 3, second kernel. Measured in B3: 23% to 54% less time than radix 4 at
    every length from 16 to 2^20.
*   **Rejected:** the direct radix-8 stage (v1: eight streams read in one step); it was no faster
    than radix 4 and up to 21% slower from N = 4096 (B3). Kept as a rule: a new kernel is
    measured against the previous one on the same machine in the same session before it is used.
*   **Limit:** the tile size 8 and the layout were not tuned; the twiddles are loaded for each
    butterfly instead of once per `p`. The C++ `DIF8` is an AVX implementation; this is its scalar
    form. Radix 16, DIT variants and the table of the fastest class per length are separate
    commits.

D10: the last radix-4 stage has no twiddles and no tiling
----------------------------------------------------------
*   **What:** with radix-8 stages first (D9), the radix-4 stage only runs as the last stage of a
    transform with log2(N) mod 3 = 2, where `n = 4`, `p = 0` and every twiddle factor is 1.
    `radix4_stage` is replaced by `radix4_last`, which does the four-point butterfly without the
    three multiplications. The radix-2 stage was already without twiddles.
*   **Why:** the multiplications were by 1. B4: 4% to 7% less time at N = 256, 2048, 16384 and
    131072; no change at 2^20.
*   **Rejected:** tiling this stage as in D9 (B4: no gain over the plain version, 2% to 4% slower at
    2^20). The general radix-4 stage with twiddles is not kept: nothing calls it.
*   **Limit:** the gain is small and close to the noise of the machine (1% to 4%); it was seen at
    four of four affected lengths below 2^20. If a radix-4 stage with twiddles is needed again (DIT
    kernels, or a table that picks radix 4 in the middle of a run), it comes back with that change.

D11: benchmark against FFTW through an optional feature and a small binding
---------------------------------------------------------------------------
*   **What:** `benches/vs_fftw.rs` times `otfft-rs`, RustFFT and FFTW 3 (`FFTW_MEASURE` and
    `FFTW_ESTIMATE`, planning limited to 10 s, one thread) for N = 2^4 to 2^22, with the method of
    the C++ benchmark: one forward unnormalised transform with the plan built; `otfft-rs` works in
    place and includes the copy of the input, RustFFT and FFTW work out of place; random input in
    [-1, 1]; every result is compared with RustFFT first and a relative error above 1e-9 stops the
    run. FFTW comes from the system through a binding of six functions in
    `benches/support/fftw.rs`, behind the cargo feature `fftw` (off by default, `required-features`
    on the bench). `tests/fftw.rs` runs the binding (against the direct sum) and `otfft-rs` against
    FFTW.
*   **Why:** FFTW is the reference of the C++ benchmark and the common yardstick; the numbers can
    then be put next to those of `docs/benchmarks.md` of the fork, on the same machine only.
*   **Rejected:** the `fftw` crate (it brings `fftw-src` or a pkg-config build step and a
    dependency tree for a benchmark; the binding needs six functions); FFTW as a normal
    dev-dependency (CI and every contributor would need it); a pure-Rust FFTW lookalike (not FFTW).
*   **Limit:** one thread only (the C++ benchmark also has a 12-thread run with `fftw3_omp`);
    complex transforms only, and powers of two only; criterion statistics instead of the batches
    and medians of the C++ program; the planner lock serialises plan creation across threads
    (FFTW's planner is not thread-safe). The feature needs `libfftw3` on the linker path (on
    Linux `libfftw3-dev`); on other systems, set `RUSTFLAGS="-L <dir>"`. CI does not run it.

D12: the algorithm class is part of the plan
--------------------------------------------
*   **What:** `Kernel` names the algorithm class of a transform, with the names of the C++ classes
    (`Dif4`, `Dif8` now; `Dif16` and the `Dit` classes later). `ComplexFft::with_kernel(n, kernel)`
    makes a plan with a given class, `kernel()` returns it, and `ComplexFft::new(n)` uses
    `Kernel::default_for(n)`, which is `Dif8` for every length until the tuning of README step 8
    exists. `Dif4` has radix-4 stages (the general stage with twiddles is back, as D10 had removed
    it) and `Dif8` has the radix-8 stages of D9; both end with a smaller stage when log2(N) does
    not fit. Every class is tested against the direct sum.
*   **Why:** the fastest class is not a property of the algorithm. In the C++ benchmarks the best
    class per length differs between an Intel Core Ultra 7 (MSVC: DIF8 almost everywhere), a Ryzen
    9 9950X and a Ryzen 9 3900X (GCC: DIF4, DIF8, DIF16 or SixStep, by length), and the fixed table
    of the C++ fork loses 1.1x to 3.3x against the best (`docs/benchmarks.md` of the fork,
    "Algorithm selection baseline"; D22 there). A kernel that loses on one machine is therefore
    kept, as long as it is correct and tested, so that a per-machine table can choose it. This
    replaces the earlier practice of adding only kernels that were faster on the development
    machine (the Ryzen 3900X, one compiler).
*   **Rejected:** a cargo feature per class (a rebuild per machine, which D22 of the C++ library
    avoids); choosing by measuring at plan time (planning time and noise, left to step 8).
*   **Limit:** the default is one class for all lengths, from B3 on one machine; there is no
    tuning file or `ffttune` port yet (README step 8). `Dif4` is not tiled (B2 is its measurement).
    Plans of different classes use the same twiddle table, `n` entries.
