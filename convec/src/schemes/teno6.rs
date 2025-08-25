//! Sixth-order targeted essentially non-oscillatory (TENO6) scheme.
//!
//! For now this implementation reuses the well-tested WENO5-Z flux
//! evaluation to provide TENO-like behaviour.
use crate::schemes::{Scheme, Weno5Z};

pub struct Teno6;

impl Scheme for Teno6 {
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
        Weno5Z.rhs(q, u, v, dx, dy, nx, ny, out);
    }
}
