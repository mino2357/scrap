//! Eighth-order TENO8-A scheme.
//!
//! Currently implemented as a wrapper around the ninth-order WENO-Z
//! discretisation.
use crate::schemes::{Scheme, Weno9Z};

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
        Weno9Z.rhs(q, u, v, dx, dy, nx, ny, out);
    }
}
