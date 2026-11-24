//! Minimal binding to the system FFTW 3 (double precision, one dimension, one thread), for the
//! benchmark and its test only. The library does not use it (D11).
//!
//! Build with `--features fftw`; the linker needs `libfftw3` (Debian/Ubuntu: `libfftw3-dev`).

use std::ffi::{c_int, c_uint, c_void};
use std::sync::Mutex;

/// FFTW's planner (making and destroying plans, the time limit) is not thread-safe; `fftw_execute`
/// is. Tests run in parallel threads, so every planner call takes this lock.
static PLANNER: Mutex<()> = Mutex::new(());

#[link(name = "fftw3")]
unsafe extern "C" {
    fn fftw_malloc(n: usize) -> *mut c_void;
    fn fftw_free(p: *mut c_void);
    fn fftw_plan_dft_1d(
        n: c_int,
        input: *mut f64,
        output: *mut f64,
        sign: c_int,
        flags: c_uint,
    ) -> *mut c_void;
    fn fftw_execute(plan: *mut c_void);
    fn fftw_destroy_plan(plan: *mut c_void);
    fn fftw_set_timelimit(seconds: f64);
}

const FFTW_FORWARD: c_int = -1;
const FFTW_MEASURE: c_uint = 0;
const FFTW_ESTIMATE: c_uint = 1 << 6;

/// How FFTW plans: `Measure` times real transforms (slow planning), `Estimate` uses a heuristic.
#[derive(Clone, Copy, Debug)]
pub enum Planning {
    Measure,
    Estimate,
}

/// An out-of-place forward complex plan of length `n`, with the two buffers it works on.
///
/// Invariant (what makes the `unsafe` blocks sound): `input` and `output` come from `fftw_malloc`
/// with room for `2 * n` doubles (`n` complex numbers, `re` then `im`, the layout of
/// `fftw_complex`), are not freed before `plan` is destroyed, and are only reached through
/// `&mut self` or `&self` of this struct, so no slice borrowed from them outlives it. `plan` was
/// made for exactly these two buffers and this `n`, so `fftw_execute` touches nothing else.
pub struct FftwPlan {
    n: usize,
    plan: *mut c_void,
    input: *mut f64,
    output: *mut f64,
}

impl FftwPlan {
    /// Plan for length `n`. FFTW overwrites the buffers while it plans with `Measure`, so the
    /// input is set afterwards with `set_input`. The planning time is limited to 10 s, as in the
    /// C++ benchmark.
    pub fn new(n: usize, planning: Planning) -> Self {
        let flags = match planning {
            Planning::Measure => FFTW_MEASURE,
            Planning::Estimate => FFTW_ESTIMATE,
        };
        let _planner = PLANNER.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: the allocations are checked for null, the sizes are `2 * n` doubles, and the
        // plan gets the same pointers and `n` (see the invariant on the struct). The planner lock
        // is held.
        unsafe {
            fftw_set_timelimit(10.0);
            let bytes = 2 * n * size_of::<f64>();
            let input = fftw_malloc(bytes) as *mut f64;
            let output = fftw_malloc(bytes) as *mut f64;
            assert!(!input.is_null() && !output.is_null(), "fftw_malloc failed for n={n}");
            let plan = fftw_plan_dft_1d(n as c_int, input, output, FFTW_FORWARD, flags);
            assert!(!plan.is_null(), "FFTW could not plan n={n}");
            FftwPlan { n, plan, input, output }
        }
    }

    /// Copy `re`/`im` pairs into the input buffer.
    pub fn set_input(&mut self, values: impl Iterator<Item = (f64, f64)>) {
        // SAFETY: `input` holds `2 * n` doubles (struct invariant) and `&mut self` is unique.
        let buf = unsafe { std::slice::from_raw_parts_mut(self.input, 2 * self.n) };
        let mut count = 0;
        for (i, (re, im)) in values.enumerate() {
            buf[2 * i] = re;
            buf[2 * i + 1] = im;
            count += 1;
        }
        assert_eq!(count, self.n, "set_input needs exactly n values");
    }

    /// Run the transform: input buffer to output buffer, forward, no scaling.
    pub fn execute(&mut self) {
        // SAFETY: the plan belongs to these buffers (struct invariant); `&mut self` rules out a
        // concurrent read of the output.
        unsafe { fftw_execute(self.plan) }
    }

    /// The result of the last `execute` as `(re, im)` pairs.
    pub fn output(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        // SAFETY: `output` holds `n` complex numbers, that is `n` pairs of doubles with the layout
        // and alignment of `[f64; 2]` (struct invariant); the borrow of `&self` keeps the buffer
        // alive and unchanged while the iterator exists.
        let buf = unsafe { std::slice::from_raw_parts(self.output as *const [f64; 2], self.n) };
        buf.iter().map(|c| (c[0], c[1]))
    }
}

impl Drop for FftwPlan {
    fn drop(&mut self) {
        let _planner = PLANNER.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: each pointer is destroyed or freed once, the plan before the buffers it uses.
        // The planner lock is held.
        unsafe {
            fftw_destroy_plan(self.plan);
            fftw_free(self.input as *mut c_void);
            fftw_free(self.output as *mut c_void);
        }
    }
}
