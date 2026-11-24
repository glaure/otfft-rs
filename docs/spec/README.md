Behaviour spec
==============

What `otfft-rs` should do, in our own words. It describes what a caller sees and gets from the
public API in `src/lib.rs`, not how the code is organised. `otfft-rs` is a port of OTFFT by OK
Ojisan (Takuya OKAHISA, <http://wwwa.pikara.ne.jp/okojisan/otfft-en/download.html>) through the
DEWETRON fork (<https://github.com/DEWETRON/otfft>). The C++ spec of the fork (`docs/spec/` there)
is the reference for behaviour; the rules here say where this port is the same and where it differs.

How to read it
--------------

Each rule has an id (`AREA-n`) and a status tag:

*   `[done]`: built, and covered by a test. The test is named (`tests/dft.rs`).
*   `[built]`: built, but known only from reading the code. No test checks it.
*   `[gap]`: the code does something else or something poorer. The rule says what is right.
*   `[planned]`: decided, not built.

Each area file ends with three lists: **In the C++ library** (what we keep or change from it),
**Decided** (links to `docs/decisions.md`) and **Questions**.

How to change it
----------------

*   Describe what the caller sees and gets, not how the code works.
*   A rule about the C++ library goes in only after someone has read or run it. Do not guess.
*   When a rule's status changes, change it in the same commit as the code or test.
*   Behaviour that differs from the C++ library is recorded in `docs/decisions.md`.

Files
-----

| File         | Area                                                              |
| ------------ | ----------------------------------------------------------------- |
| `complex.md` | Complex in-place FFT: plan, forward, inverse, normalisation       |
| `errors.md`  | Refused lengths and buffers                                       |
| `build.md`   | Crate, dependencies, tests, platforms                             |
| `planned.md` | Areas of the C++ library that are not ported yet                  |
