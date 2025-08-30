//! Eighth-order TENO8-A scheme.
//!
//! Currently an alias of the TENO9-A implementation: both use the
//! WENO9 family candidates with the same TENO cut-off. Kept to retain
//! configuration compatibility; may diverge when an 8th-order optimal
//! coefficient set is added.
use crate::schemes::{Scheme, Teno9A};

pub struct Teno8A;

impl Scheme for Teno8A {
    fn rhs(
        &self,
        q: &[f64],
        u: &[f64],
        v: &[f64],
        dx: f64,
        dy: f64,
        nx: usize,
        ny: usize,
        out: &mut [f64],
    ) {
        // Delegate to TENO9-A. Behaviour is identical in the current codebase.
        Teno9A.rhs(q, u, v, dx, dy, nx, ny, out);
    }
}
