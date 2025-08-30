//! Ninth-order TENO9-A scheme.
//!
//! Implemented as a TENO cut-off on the WENO9 family reconstruction.
//! In smooth regions it reduces towards the optimal linear combination; near
//! discontinuities it discards non-smooth substencils and renormalises the
//! remaining ones. 数値流束は面速度の符号に基づく upwind 面流束。
use crate::schemes::Scheme;
use crate::utils::{idx, pid};
use super::weno9_core::{cvals_betas, D, EPS};

pub struct Teno9A;


fn teno9a_reconstruct(arr: &[f64; 9]) -> f64 {
    let (cval, beta) = cvals_betas(arr);
    let tau = (beta[0] - beta[4]).abs();

    // WENO-Z baseline
    let mut alpha = [0.0; 5];
    let mut sum = 0.0;
    for k in 0..5 {
        alpha[k] = D[k] * (1.0 + (tau / (beta[k] + EPS)).powi(2));
        sum += alpha[k];
    }
    let mut wz = [0.0; 5];
    for k in 0..5 {
        wz[k] = alpha[k] / sum;
    }

    // TENO cut-off
    let ct = 1e-7;
    let q = 6;
    let mut a = [0.0; 5];
    let mut s = 0.0;
    for k in 0..5 {
        let rk = ((tau + EPS) / (beta[k] + EPS)).powi(q);
        a[k] = if rk < ct { wz[k] } else { 0.0 };
        s += a[k];
    }
    if s > 0.0 {
        a.iter_mut().for_each(|x| *x /= s);
        a.iter().zip(cval.iter()).map(|(w, c)| w * c).sum()
    } else {
        wz.iter().zip(cval.iter()).map(|(w, c)| w * c).sum()
    }
}

fn teno9a_fd_flux_faces_1d_upwind(q: &[f64], vel: &[f64]) -> Vec<f64> {
    let n = q.len();
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let ip = pid(i as isize + 1, n);
        let u = 0.5 * (vel[i] + vel[ip]);
        let mut arrp = [0.0; 9];
        for k in 0..9 { arrp[k] = q[pid(i as isize + k as isize - 4, n)]; }
        let ql = teno9a_reconstruct(&arrp);
        let mut arrm = [0.0; 9];
        for k in 0..9 { arrm[k] = q[pid(i as isize + 5 - k as isize, n)]; }
        let qr = teno9a_reconstruct(&arrm);
        let qu = if u >= 0.0 { ql } else { qr };
        fh[i] = u * qu;
    }
    fh
}

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
        let mut dfx = vec![0.0; nx * ny];
        for j in 0..ny {
            let mut qq = vec![0.0; nx];
            let mut uu = vec![0.0; nx];
            for i in 0..nx {
                let k = idx(i, j, nx);
                qq[i] = q[k];
                uu[i] = u[k];
            }
            let fh = teno9a_fd_flux_faces_1d_upwind(&qq, &uu);
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
            let gh = teno9a_fd_flux_faces_1d_upwind(&qq, &vv);
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
