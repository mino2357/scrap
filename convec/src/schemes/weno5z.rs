//! 5次精度の加重基本的非振動 (WENO) スキームを用いた移流項評価。
//! 実装は Borges らによる改良型の WENO-Z スキーム\[1\]に基づいている。
//!
//! [1] R. Borges, M. Carmona, B. Costa, and W. S. Don,
//!     "An improved weighted essentially non-oscillatory scheme for
//!     hyperbolic conservation laws", *Journal of Computational Physics*,
//!     227(6), 3191-3211, 2008.
use crate::schemes::Scheme;
use crate::utils::{idx, pid};
use super::weno5_core::{cvals_betas_5, D5, EPS5};

pub struct Weno5Z;

fn weno5z_fd_flux_faces_1d(f: &[f64], q: &[f64], alpha: f64) -> Vec<f64> {
    let n = f.len();
    let mut fp = vec![0.0; n];
    let mut fm = vec![0.0; n];
    for i in 0..n {
        fp[i] = (1.0 / 2.0) * (f[i] + alpha * q[i]);
        fm[i] = (1.0 / 2.0) * (f[i] - alpha * q[i]);
    }
    let mut fh = vec![0.0; n];
    for i in 0..n {
        // positive flux
        let arrp = [
            fp[pid(i as isize - 2, n)],
            fp[pid(i as isize - 1, n)],
            fp[i],
            fp[pid(i as isize + 1, n)],
            fp[pid(i as isize + 2, n)],
        ];
        let (cvalp, betap) = cvals_betas_5(&arrp);
        let taup = (betap[0] - betap[2]).abs();
        let a0 = D5[0] * (1.0 + (taup / (EPS5 + betap[0])).powi(2));
        let a1 = D5[1] * (1.0 + (taup / (EPS5 + betap[1])).powi(2));
        let a2 = D5[2] * (1.0 + (taup / (EPS5 + betap[2])).powi(2));
        let sp = a0 + a1 + a2;
        let w0 = a0 / sp;
        let w1 = a1 / sp;
        let w2 = a2 / sp;
        let fph = w0 * cvalp[0] + w1 * cvalp[1] + w2 * cvalp[2];

        // negative flux (mirror stencil about i+1/2)
        let arrm = [
            fm[pid(i as isize + 3, n)],
            fm[pid(i as isize + 2, n)],
            fm[pid(i as isize + 1, n)],
            fm[i],
            fm[pid(i as isize - 1, n)],
        ];
        let (cvalm, betam) = cvals_betas_5(&arrm);
        let taum = (betam[0] - betam[2]).abs();
        let a0m = D5[0] * (1.0 + (taum / (EPS5 + betam[0])).powi(2));
        let a1m = D5[1] * (1.0 + (taum / (EPS5 + betam[1])).powi(2));
        let a2m = D5[2] * (1.0 + (taum / (EPS5 + betam[2])).powi(2));
        let sm = a0m + a1m + a2m;
        let w0m = a0m / sm;
        let w1m = a1m / sm;
        let w2m = a2m / sm;
        let fmh = w2m * cvalm[0] + w1m * cvalm[1] + w0m * cvalm[2];
        fh[i] = fph + fmh;
    }
    fh
}

impl Scheme for Weno5Z {
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
            let mut f = vec![0.0; nx];
            let mut qq = vec![0.0; nx];
            let mut amax = 0.0;
            for i in 0..nx {
                let k = idx(i, j, nx);
                f[i] = u[k] * q[k];
                qq[i] = q[k];
                let a = u[k].abs();
                if a > amax {
                    amax = a;
                }
            }
            let fh = weno5z_fd_flux_faces_1d(&f, &qq, amax);
            for i in 0..nx {
                let im = pid(i as isize - 1, nx);
                dfx[idx(i, j, nx)] = (fh[i] - fh[im]) / dx;
            }
        }
        let mut dfy = vec![0.0; nx * ny];
        for i in 0..nx {
            let mut g = vec![0.0; ny];
            let mut qq = vec![0.0; ny];
            let mut amax = 0.0;
            for j in 0..ny {
                let k = idx(i, j, nx);
                g[j] = v[k] * q[k];
                qq[j] = q[k];
                let a = v[k].abs();
                if a > amax {
                    amax = a;
                }
            }
            let gh = weno5z_fd_flux_faces_1d(&g, &qq, amax);
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
