use crate::schemes::Scheme;
use crate::utils::{idx, pid};

/// 一次精度の風上差分スキーム。
pub struct Upwind1;

impl Scheme for Upwind1 {
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
        for j in 0..ny {
            for i in 0..nx {
                let k = idx(i, j, nx);
                let im = pid(i as isize - 1, nx);
                let ip = pid(i as isize + 1, nx);
                let jm = pid(j as isize - 1, ny);
                let jp = pid(j as isize + 1, ny);
                let dqx = if u[k] >= 0.0 {
                    (q[k] - q[idx(im, j, nx)]) / dx
                } else {
                    (q[idx(ip, j, nx)] - q[k]) / dx
                };
                let dqy = if v[k] >= 0.0 {
                    (q[k] - q[idx(i, jm, nx)]) / dy
                } else {
                    (q[idx(i, jp, nx)] - q[k]) / dy
                };
                out[k] = -(u[k] * dqx + v[k] * dqy);
            }
        }
    }
}
