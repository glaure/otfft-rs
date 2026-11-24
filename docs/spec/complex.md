Complex FFT
===========

`ComplexFft`: in-place discrete Fourier transform on `&mut [Complex]` (`Complex`: two `f64`, `re`
and `im`). Power-of-two lengths only for now (see `planned.md`).

Rules
-----

*   **CFFT-1** `[done]` `ComplexFft::new(n)` makes a plan for length `n`, a power of two from 1
    up. (`forward_matches_the_direct_sum`, lengths 1 to 2^11)
*   **CFFT-2** `[done]` A forward transform replaces `x` by its spectrum
    `X[k] = sum over j of x[j] * exp(-2 pi i j k / N)`, before scaling. It matches the direct sum
    with an error below `1e-9 * N`. (`forward_matches_the_direct_sum`)
*   **CFFT-3** `[done]` An inverse transform is the same sum with `exp(+...)`, before scaling.
    (`inverse_matches_the_direct_sum`, lengths 1 to 2^10)
*   **CFFT-4** `[done]` Forward followed by inverse gives back the input when the scales of the two
    fit: `fwd` then `inv`, or `fwd0` then `invn`. (`fwd_then_inv_returns_the_signal`, lengths 1 to
    2^16)
*   **CFFT-5** `[done]` Scales are those of the C++ library (its CFFT-4):

    | Method          | Direction | Scale     |
    | --------------- | --------- | --------- |
    | `fwd_alloc`     | forward   | 1/N       |
    | `fwd0_alloc`    | forward   | none      |
    | `fwdu_alloc`    | forward   | 1/sqrt(N) |
    | `inv_alloc`     | inverse   | none      |
    | `invu_alloc`    | inverse   | 1/sqrt(N) |
    | `invn_alloc`    | inverse   | 1/N       |

    `process(x, scratch, Direction, Scale)` takes the direction and the scale as arguments; the
    methods above are shortcuts for it. (`the_scales_are_those_of_the_cpp_library`)
*   **CFFT-6** `[done]` A pure tone with `f` cycles over the length gives a spectrum with the value 1
    in bin `f` and 0 elsewhere after `fwd`. (`a_pure_tone_lands_in_one_bin`)
*   **CFFT-7** `[done]` A plan holds only read-only data (the twiddle factors), is `Send + Sync` and
    can be used from several threads at the same time, each with its own buffers. This differs from
    the C++ library (its FAC-6). (`a_plan_can_be_shared_between_threads`)
*   **CFFT-8** `[done]` `process` needs a scratch buffer of `scratch_len()` elements (equal to `N`).
    The content of the scratch buffer on entry and on return is not defined, and one buffer can be
    reused for any number of calls. The `*_alloc` methods allocate one per call.
    (`one_scratch_buffer_serves_many_calls`)
*   **CFFT-9** `[done]` `Complex` has the memory layout of `OTFFT::complex_t` of the C++ library:
    16 bytes, `re` first. A buffer can be handed to C++ code without copying.
    (`complex_has_the_layout_of_the_cpp_complex_t`)
*   **CFFT-10** `[built]` No `unsafe` code in the transform, and no allocation in `process`.
*   **CFFT-11** `[gap]` Speed: scalar radix 2, much slower than the C++ library and than RustFFT.
    Right: at least as fast as the C++ library at the same length (planned steps 2 to 4 in the
    README). Not measured yet.

In the C++ library
------------------

*   Same names and scales: `fwd` divides by N, `inv` does not (kept, D2).
*   The C++ object owns the work buffer and must not be shared; here the caller brings it (D3).
*   The C++ `fwd`/`inv` take `complex_vector` of the right length without checking; here the length
    is checked (`errors.md`).

Decided
-------

*   D2, D3, D5.

Questions
---------

*   Should `process` also accept a scratch buffer longer than `N`? Today it must be exactly `N`.
