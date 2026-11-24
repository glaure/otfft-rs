Rules for agents in this repository
===================================

The key words MUST, MUST NOT, SHOULD, SHOULD NOT and MAY are used as in
RFC 2119 and RFC 8174. When you depart from a SHOULD, say so and why in the
commit message.

Changing code
-------------

1.  You MUST NOT change code in `src/` without a test in `tests/` (or a unit
    test in the same file) that covers the change. Transforms are checked
    against the direct DFT sum, not only against another fast transform.
2.  You SHOULD keep changes small: one step of the plan in `README.md` per
    commit, and for a kernel ported from the C++ library, one kernel per
    commit.
3.  You MUST keep the library free of dependencies. Benchmark and test helpers
    (`criterion`, `rustfft`, `realfft`) go in `[dev-dependencies]` only.
4.  You MUST NOT use `unsafe` without a comment that states the invariant, and
    a test that runs the unsafe path. Prefer safe code and let the optimiser
    work first; SIMD with `std::arch` is `unsafe` by nature (plan step 4) and
    needs the run-time feature check next to it.
5.  You MUST keep the public API names and scales of the C++ library (`fwd`
    divides by N, `inv` does not) or record the break in `docs/decisions.md`.
6.  You MUST NOT add a dependency or a `publish` setting that makes the crate
    publishable: it is private (`publish = false`, D4).

Docs
----

7.  When behaviour changes, you MUST update `docs/spec/` in the same commit,
    including the status tag and the name of the test that covers the rule.
8.  A new design choice MUST get a new `Dn` entry at the end of
    `docs/decisions.md`. Replaced entries are marked superseded, not deleted.
9.  You MUST NOT claim a speed in a doc or a commit message without a
    measurement. Measurements go in `docs/benchmarks.md` with the machine,
    the compiler version and the command; compare with the numbers of the C++
    library in the DEWETRON repository only for the same machine.

Tests and commits
-----------------

10. You MUST run `cargo test` and `cargo clippy --all-targets` before you say
    work is done. Both must report nothing to fix.
11. You MUST NOT commit or push unless asked.
12. You MUST NOT add a `Co-Authored-By:` line, or any other agent attribution
    (such as "Generated with Claude Code"), to commit messages.

This repository
---------------

`otfft-rs`: a Rust port of OTFFT (OK Ojisan), through the DEWETRON fork of the
C++ library: <https://github.com/DEWETRON/otfft>. Private project. MIT license;
the notices of the original author and of the fork stay in `LICENSE`.

### Layout

*   `src/lib.rs`: the public API and the transforms.
*   `tests/`: integration tests (`dft.rs`: transforms against the direct sum).
*   `docs/spec/`: behaviour spec with rule ids and status tags.
    `docs/decisions.md`: why it is built this way.
*   `README.md`: status and the plan, step by step.

### Checks

```bash
cargo test
cargo clippy --all-targets
```

### Reference

The C++ library is the reference for behaviour and for speed. Its spec
(`docs/spec/`), decisions (`docs/decisions.md`) and measurements
(`docs/benchmarks.md`) are in the DEWETRON repository. A rule about the C++
library goes into this spec only after someone has read or run it.
