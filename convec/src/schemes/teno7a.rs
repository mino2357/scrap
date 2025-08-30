//! Seventh-order Targeted ENO (TENO7-A) scheme.
//!
//! This implementation follows the adaptive-order formulation where
//! candidate substencils are either included with their optimal linear
//! weight or completely discarded based on a smoothness measurement.
//! The reconstruction reuses the seventh-order optimal polynomials from
//! the WENO7 family and applies the TENO cut-off (Fu et al., 2016).
//! 数値流束は面速度の符号で左右再構成から上流値を選ぶ upwind 面流束。

use crate::schemes::Scheme;
use crate::utils::{idx, pid};

pub struct Teno7A;

fn beta(a1: f64, a2: f64, a3: f64) -> f64 {
    // Smoothness indicator used by both WENO7-Z and the TENO7-A scheme.
    a1 * a1 + 0.5 * a1 * a3 + (13.0 / 3.0) * a2 * a2 + (3129.0 / 80.0) * a3 * a3
}

fn reconstruct(arr: &[f64; 7]) -> f64 {
    let im3 = arr[0];
    let im2 = arr[1];
    let im1 = arr[2];
    let i0 = arr[3];
    let ip1 = arr[4];
    let ip2 = arr[5];
    let ip3 = arr[6];

    // Linear reconstructions on each candidate stencil.
    let c0 = (-5.0 * im3 + 21.0 * im2 - 35.0 * im1 + 35.0 * i0) / 16.0;
    let c1 = (im2 - 5.0 * im1 + 15.0 * i0 + 5.0 * ip1) / 16.0;
    let c2 = (-im1 + 9.0 * i0 + 9.0 * ip1 - ip2) / 16.0;
    let c3 = (5.0 * i0 + 15.0 * ip1 - 5.0 * ip2 + ip3) / 16.0;

    // Smoothness indicators for each substencil.
    let a1 = (-1.0 / 3.0) * im3 + (3.0 / 2.0) * im2 - 3.0 * im1 + 11.0 / 6.0 * i0;
    let a2 = (-1.0 / 2.0) * im3 + 2.0 * im2 - (5.0 / 2.0) * im1 + i0;
    let a3 = (-1.0 / 6.0) * im3 + (1.0 / 2.0) * im2 - (1.0 / 2.0) * im1 + 1.0 / 6.0 * i0;
    let b0 = beta(a1, a2, a3);
    let a1 = (1.0 / 6.0) * im2 - im1 + (1.0 / 2.0) * i0 + 1.0 / 3.0 * ip1;
    let a2 = (1.0 / 2.0) * im1 - i0 + (1.0 / 2.0) * ip1;
    let a3 = (-1.0 / 6.0) * im2 + (1.0 / 2.0) * im1 - (1.0 / 2.0) * i0 + 1.0 / 6.0 * ip1;
    let b1 = beta(a1, a2, a3);
    let a1 = (-1.0 / 3.0) * im1 - (1.0 / 2.0) * i0 + ip1 - 1.0 / 6.0 * ip2;
    let a2 = (1.0 / 2.0) * im1 - i0 + (1.0 / 2.0) * ip1;
    let a3 = (-1.0 / 6.0) * im1 + (1.0 / 2.0) * i0 - (1.0 / 2.0) * ip1 + 1.0 / 6.0 * ip2;
    let b2 = beta(a1, a2, a3);
    let a1 = (-11.0 / 6.0) * i0 + 3.0 * ip1 - (3.0 / 2.0) * ip2 + 1.0 / 3.0 * ip3;
    let a2 = i0 - (5.0 / 2.0) * ip1 + 2.0 * ip2 - (1.0 / 2.0) * ip3;
    let a3 = (-1.0 / 6.0) * i0 + (1.0 / 2.0) * ip1 - (1.0 / 2.0) * ip2 + 1.0 / 6.0 * ip3;
    let b3 = beta(a1, a2, a3);

    // Global smoothness measure and linear weights.
    let eps = 1e-6;
    let tau7 = (b0 - b3).abs();
    let d0 = 1.0 / 35.0;
    let d1 = 12.0 / 35.0;
    let d2 = 18.0 / 35.0;
    let d3 = 4.0 / 35.0;

    // TENO cut-off parameters.
    let ct = 1e-6;
    let q = 6; // power recommended in the literature

    // Nonlinear weights following the WENO7-Z formulation.
    let alpha0 = d0 * (1.0 + (tau7 / (b0 + eps)).powi(q));
    let alpha1 = d1 * (1.0 + (tau7 / (b1 + eps)).powi(q));
    let alpha2 = d2 * (1.0 + (tau7 / (b2 + eps)).powi(q));
    let alpha3 = d3 * (1.0 + (tau7 / (b3 + eps)).powi(q));
    let alpha_sum = alpha0 + alpha1 + alpha2 + alpha3;
    let w0z = alpha0 / alpha_sum;
    let w1z = alpha1 / alpha_sum;
    let w2z = alpha2 / alpha_sum;
    let w3z = alpha3 / alpha_sum;

    // Smoothness sensors for each substencil.  The ratio places the global
    // smoothness measure in the numerator so that the sensor is small in smooth
    // regions (tau ≪ beta) and O(1) near discontinuities (tau ≈ beta).
    let r0 = ((tau7 + eps) / (b0 + eps)).powi(q);
    let r1 = ((tau7 + eps) / (b1 + eps)).powi(q);
    let r2 = ((tau7 + eps) / (b2 + eps)).powi(q);
    let r3 = ((tau7 + eps) / (b3 + eps)).powi(q);

    // Include or discard substencils according to the sensor.
    let a0 = if r0 < ct { w0z } else { 0.0 };
    let a1 = if r1 < ct { w1z } else { 0.0 };
    let a2 = if r2 < ct { w2z } else { 0.0 };
    let a3 = if r3 < ct { w3z } else { 0.0 };
    let s = a0 + a1 + a2 + a3;
    let (w0, w1, w2, w3) = if s > 0.0 {
        (a0 / s, a1 / s, a2 / s, a3 / s)
    } else {
        // All stencils flagged as non-smooth; revert to WENO-Z weights.
        (w0z, w1z, w2z, w3z)
    };

    w0 * c0 + w1 * c1 + w2 * c2 + w3 * c3
}

fn teno7a_fd_flux_faces_1d_upwind(q: &[f64], vel: &[f64]) -> Vec<f64> {
    let n = q.len();
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let ip = pid(i as isize + 1, n);
        let u = 0.5 * (vel[i] + vel[ip]);
        let arrp = [
            q[pid(i as isize - 3, n)],
            q[pid(i as isize - 2, n)],
            q[pid(i as isize - 1, n)],
            q[i],
            q[ip],
            q[pid(i as isize + 2, n)],
            q[pid(i as isize + 3, n)],
        ];
        let ql = reconstruct(&arrp);
        let arrm = [
            q[pid(i as isize + 4, n)],
            q[pid(i as isize + 3, n)],
            q[pid(i as isize + 2, n)],
            q[ip],
            q[i],
            q[pid(i as isize - 1, n)],
            q[pid(i as isize - 2, n)],
        ];
        let qr = reconstruct(&arrm);
        let qu = if u >= 0.0 { ql } else { qr };
        fh[i] = u * qu;
    }
    fh
}

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
        let mut dfx = vec![0.0; nx * ny];
        for j in 0..ny {
            let mut qq = vec![0.0; nx];
            let mut uu = vec![0.0; nx];
            for i in 0..nx {
                let k = idx(i, j, nx);
                qq[i] = q[k];
                uu[i] = u[k];
            }
            let fh = teno7a_fd_flux_faces_1d_upwind(&qq, &uu);
            for i in 0..nx {
                let im = pid(i as isize - 1, nx);
                dfx[idx(i, j, nx)] = (fh[i] - fh[im]) / dx;
            }
        }
        let mut dfy = vec![0.0; nx * ny];
        for i in 0..nx {
            let mut qq = vec![0.0; ny];
            let mut vv = vec![0.0; ny];
            for j in 0..ny {
                let k = idx(i, j, nx);
                qq[j] = q[k];
                vv[j] = v[k];
            }
            let gh = teno7a_fd_flux_faces_1d_upwind(&qq, &vv);
            for j in 0..ny {
                let jm = pid(j as isize - 1, ny);
                dfy[idx(i, j, nx)] = (gh[j] - gh[jm]) / dy;
            }
        }
        for k in 0..nx * ny {
            out[k] = -(dfx[k] + dfy[k]);
        }
    }
}
