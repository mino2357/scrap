//! Ninth-order TENO9-A scheme.
//!
//! Uses the existing WENO9-Z implementation for flux evaluation.
use crate::schemes::{Scheme, Weno9Z};

pub struct Teno9A;

impl Scheme for Teno9A {
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
        Weno9Z.rhs(q, u, v, dx, dy, nx, ny, out);
    }
}
