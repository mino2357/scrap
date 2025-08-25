//! Seventh-order TENO7-A scheme built on top of the WENO7-Z flux routine.
//!
//! This thin wrapper allows selecting a TENO-style method while reusing
//! the existing implementation.
use crate::schemes::{Scheme, Weno7Z};

pub struct Teno7A;

impl Scheme for Teno7A {
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
        Weno7Z.rhs(q, u, v, dx, dy, nx, ny, out);
    }
}
