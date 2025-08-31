//! 14次精度の中心差分による移流項の離散化を行うスキーム。
//! 係数は Fornberg による有限差分公式の一般化[1]に基づく。
//!
//! [1] B. Fornberg, "Generation of finite difference formulas on arbitrarily spaced grids",
//! *Mathematics of Computation*, 51(184), 699-706, 1988.
use crate::schemes::Scheme;
use crate::utils::{idx, pid};

pub struct Centered14;

// Coefficients for the 14th-order centered finite-difference approximation of the first derivative.
// Generated via SymPy's `finite_diff_weights` following Fornberg's algorithm.
const C14: [f64; 15] = [
    -1.0 / 24024.0,
    7.0 / 10296.0,
    -7.0 / 1320.0,
    7.0 / 264.0,
    -7.0 / 72.0,
    7.0 / 24.0,
    -7.0 / 8.0,
    0.0,
    7.0 / 8.0,
    -7.0 / 24.0,
    7.0 / 72.0,
    -7.0 / 264.0,
    7.0 / 1320.0,
    -7.0 / 10296.0,
    1.0 / 24024.0,
];

fn d1x(f: &[f64], nx: usize, ny: usize, dx: f64, out: &mut [f64]) {
    for j in 0..ny {
        for i in 0..nx {
            let mut acc = 0.0;
            for s in -7..=7 {
                if s != 0 {
                    let w = C14[(s + 7) as usize];
                    let ii = pid(i as isize + s as isize, nx);
                    acc += w * f[idx(ii, j, nx)];
                }
            }
            out[idx(i, j, nx)] = acc / dx;
        }
    }
}

fn d1y(f: &[f64], nx: usize, ny: usize, dy: f64, out: &mut [f64]) {
    for j in 0..ny {
        for i in 0..nx {
            let mut acc = 0.0;
            for s in -7..=7 {
                if s != 0 {
                    let w = C14[(s + 7) as usize];
                    let jj = pid(j as isize + s as isize, ny);
                    acc += w * f[idx(i, jj, nx)];
                }
            }
            out[idx(i, j, nx)] = acc / dy;
        }
    }
}

impl Scheme for Centered14 {
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
        let mut dqx = vec![0.0; nx * ny];
        let mut dqy = vec![0.0; nx * ny];
        d1x(q, nx, ny, dx, &mut dqx);
        d1y(q, nx, ny, dy, &mut dqy);

        let mut uq = vec![0.0; nx * ny];
        let mut vq = vec![0.0; nx * ny];
        for k in 0..nx * ny {
            uq[k] = u[k] * q[k];
            vq[k] = v[k] * q[k];
        }

        let mut dx_uq = vec![0.0; nx * ny];
        let mut dy_vq = vec![0.0; nx * ny];
        d1x(&uq, nx, ny, dx, &mut dx_uq);
        d1y(&vq, nx, ny, dy, &mut dy_vq);

        for k in 0..nx * ny {
            out[k] = -0.5 * (u[k] * dqx[k] + dx_uq[k] + v[k] * dqy[k] + dy_vq[k]);
        }
    }
}
