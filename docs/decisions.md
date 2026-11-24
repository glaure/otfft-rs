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

D5: scalar radix 2 first
------------------------
*   **What:** step 0 is a plain radix-2 Stockham with a twiddle table, power-of-two lengths only.
*   **Why:** a correct and tested baseline to measure the radix-4/8/16 and SIMD steps against, and a
    harness (direct-sum tests) that every later kernel must pass.
*   **Limit:** slow; not a replacement for the C++ library yet (CFFT-11).
