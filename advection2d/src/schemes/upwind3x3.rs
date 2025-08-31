use crate::schemes::Scheme;
use crate::utils::{idx, pid};

/// First-order upwind scheme with a 3x3 stencil.
/// The derivative is evaluated with a two-dimensional Simpson rule
/// taken in the upwind direction, using diagonal information as well.
pub struct Upwind3x3;

impl Scheme for Upwind3x3 {
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
            let jm = pid(j as isize - 1, ny);
            let jp = pid(j as isize + 1, ny);
            for i in 0..nx {
                let im = pid(i as isize - 1, nx);
                let ip = pid(i as isize + 1, nx);
                let k = idx(i, j, nx);

                let dqx = if u[k] >= 0.0 {
                    (
                        q[idx(i, jm, nx)] + 4.0 * q[idx(i, j, nx)] + q[idx(i, jp, nx)]
                            - q[idx(im, jm, nx)] - 4.0 * q[idx(im, j, nx)] - q[idx(im, jp, nx)]
                    ) / (6.0 * dx)
                } else {
                    (
                        q[idx(ip, jm, nx)] + 4.0 * q[idx(ip, j, nx)] + q[idx(ip, jp, nx)]
                            - q[idx(i, jm, nx)] - 4.0 * q[idx(i, j, nx)] - q[idx(i, jp, nx)]
                    ) / (6.0 * dx)
                };

                let dqy = if v[k] >= 0.0 {
                    (
                        q[idx(im, j, nx)] + 4.0 * q[idx(i, j, nx)] + q[idx(ip, j, nx)]
                            - q[idx(im, jm, nx)] - 4.0 * q[idx(i, jm, nx)] - q[idx(ip, jm, nx)]
                    ) / (6.0 * dy)
                } else {
                    (
                        q[idx(im, jp, nx)] + 4.0 * q[idx(i, jp, nx)] + q[idx(ip, jp, nx)]
                            - q[idx(im, j, nx)] - 4.0 * q[idx(i, j, nx)] - q[idx(ip, j, nx)]
                    ) / (6.0 * dy)
                };

                out[k] = -(u[k] * dqx + v[k] * dqy);
            }
        }
    }
}

