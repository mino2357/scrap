use crate::schemes::Scheme;
use crate::utils::{idx, pid};

/// Second-order central difference scheme with a 4x4 stencil.
/// This scheme employs the two-dimensional Simpson 3/8 rule
/// to incorporate a wider neighbourhood for improved isotropy.
pub struct Centered4x4;

impl Scheme for Centered4x4 {
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
            let jm2 = pid(j as isize - 2, ny);
            let jm1 = pid(j as isize - 1, ny);
            let jp1 = pid(j as isize + 1, ny);
            let jp2 = pid(j as isize + 2, ny);
            for i in 0..nx {
                let im2 = pid(i as isize - 2, nx);
                let im1 = pid(i as isize - 1, nx);
                let ip1 = pid(i as isize + 1, nx);
                let ip2 = pid(i as isize + 2, nx);
                let k = idx(i, j, nx);

                let dqx = (q[idx(ip1, jm1, nx)]
                    + 3.0 * q[idx(ip1, j, nx)]
                    + 3.0 * q[idx(ip1, jp1, nx)]
                    + q[idx(ip1, jp2, nx)]
                    - q[idx(im1, jm2, nx)]
                    - 3.0 * q[idx(im1, jm1, nx)]
                    - 3.0 * q[idx(im1, j, nx)]
                    - q[idx(im1, jp1, nx)])
                    / (8.0 * dx);

                let dqy = (q[idx(im1, jp1, nx)]
                    + 3.0 * q[idx(i, jp1, nx)]
                    + 3.0 * q[idx(ip1, jp1, nx)]
                    + q[idx(ip2, jp1, nx)]
                    - q[idx(im2, jm1, nx)]
                    - 3.0 * q[idx(im1, jm1, nx)]
                    - 3.0 * q[idx(i, jm1, nx)]
                    - q[idx(ip1, jm1, nx)])
                    / (8.0 * dy);

                out[k] = -(u[k] * dqx + v[k] * dqy);
            }
        }
    }
}
