Errors
======

What the crate refuses, and what the caller then sees.

Rules
-----

*   **ERR-1** `[done]` `ComplexFft::new(n)` returns `Error::UnsupportedLength(n)` for `n = 0` and for
    a length that is not a power of two (until the mixed radix and Bluestein steps exist, see
    `planned.md`). (`refuses_lengths_and_buffers_it_cannot_handle`)
*   **ERR-2** `[done]` `process` returns `Error::WrongBufferLength { expected, got }` if `x` or the
    scratch buffer does not have exactly `len()` elements. (`refuses_lengths_and_buffers_it_cannot_handle`)
*   **ERR-3** `[done]` A refused call leaves `x` unchanged. (`a_refused_call_leaves_the_data_unchanged`)
*   **ERR-4** `[built]` `Error` implements `Display` and `std::error::Error`.
*   **ERR-5** `[built]` The `*_alloc` methods allocate their own scratch buffer and so cannot hit a
    buffer length error; they panic (`expect`) only if `x` has the wrong length. That is a caller bug,
    not a runtime condition.

In the C++ library
------------------

*   A length that is not supported is not an error there: odd lengths for the real FFT and DCT are
    undocumented (C++ spec, RFFT questions). Here every unsupported case returns an error.

Decided
-------

*   D4.

Questions
---------

*   Should the `*_alloc` methods return `Result` instead of panicking on a wrong `x` length?
