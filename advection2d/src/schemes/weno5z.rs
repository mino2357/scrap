//! 5次精度の加重基本的非振動 (WENO) スキーム。
//! 実装は Borges らによる改良型の WENO-Z スキーム\[1\]に基づき，
//! 面速度の符号で上流側再構成を選ぶ upwind 面流束を用いる。
//!
//! [1] R. Borges, M. Carmona, B. Costa, and W. S. Don,
//!     "An improved weighted essentially non-oscillatory scheme for
//!     hyperbolic conservation laws", *Journal of Computational Physics*,
//!     227(6), 3191-3211, 2008.
use crate::schemes::Scheme;
use crate::utils::{idx, pid};
use super::weno5_core::{cvals_betas_5, D5, EPS5};

pub struct Weno5Z;

fn weno5z_fd_flux_faces_1d(q: &[f64], vel: &[f64]) -> Vec<f64> {
    let n = q.len();
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let ip = pid(i as isize + 1, n);
        let u = 0.5 * (vel[i] + vel[ip]);
        // Left reconstruction at i+1/2
        let arrp = [
            q[pid(i as isize - 2, n)],
            q[pid(i as isize - 1, n)],
            q[i],
            q[ip],
            q[pid(i as isize + 2, n)],
        ];
        let (cvalp, betap) = cvals_betas_5(&arrp);
        let taup = (betap[0] - betap[2]).abs();
        let a0 = D5[0] * (1.0 + (taup / (EPS5 + betap[0])).powi(2));
        let a1 = D5[1] * (1.0 + (taup / (EPS5 + betap[1])).powi(2));
        let a2 = D5[2] * (1.0 + (taup / (EPS5 + betap[2])).powi(2));
        let sp = a0 + a1 + a2;
        let ql = (a0 / sp) * cvalp[0] + (a1 / sp) * cvalp[1] + (a2 / sp) * cvalp[2];
        // Right reconstruction at i+1/2 (mirror)
        let arrm = [
            q[pid(i as isize + 3, n)],
            q[pid(i as isize + 2, n)],
            q[ip],
            q[i],
            q[pid(i as isize - 1, n)],
        ];
        let (cvalm, betam) = cvals_betas_5(&arrm);
        let taum = (betam[0] - betam[2]).abs();
        let a0m = D5[0] * (1.0 + (taum / (EPS5 + betam[0])).powi(2));
        let a1m = D5[1] * (1.0 + (taum / (EPS5 + betam[1])).powi(2));
        let a2m = D5[2] * (1.0 + (taum / (EPS5 + betam[2])).powi(2));
        let sm = a0m + a1m + a2m;
        let qr = (a2m / sm) * cvalm[0] + (a1m / sm) * cvalm[1] + (a0m / sm) * cvalm[2];
        let qu = if u >= 0.0 { ql } else { qr };
        fh[i] = u * qu;
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
            let mut qq = vec![0.0; nx];
            let mut uu = vec![0.0; nx];
            for i in 0..nx {
                let k = idx(i, j, nx);
                qq[i] = q[k];
                uu[i] = u[k];
            }
            let fh = weno5z_fd_flux_faces_1d(&qq, &uu);
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
            let gh = weno5z_fd_flux_faces_1d(&qq, &vv);
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
