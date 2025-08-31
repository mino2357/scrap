//! 5次精度の加重基本的非振動 (WENO) スキームを用いた移流項評価。
//! 実装は Jiang & Shu による古典的な WENO-JS スキーム\[1\]に基づいている。
//!
//! [1] G.-S. Jiang and C.-W. Shu, "Efficient implementation of Weighted ENO schemes",
//! *Journal of Computational Physics*, 126(1), 202-228, 1996.
use crate::schemes::Scheme;
use crate::utils::{idx, pid};
use super::weno5_core::{cvals_betas_5, D5, EPS5};

pub struct Weno5Js;

// Upwind face flux with Jiang–Shu weights (avoids LLF directional bias)
fn weno5_fd_flux_faces_1d_upwind(q: &[f64], vel: &[f64]) -> Vec<f64> {
    let n = q.len();
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let ip = pid(i as isize + 1, n);
        let u = 0.5 * (vel[i] + vel[ip]);
        // Left reconstruction at i+1/2 using WENO5-JS weights
        let arrp = [
            q[pid(i as isize - 2, n)],
            q[pid(i as isize - 1, n)],
            q[i],
            q[ip],
            q[pid(i as isize + 2, n)],
        ];
        let (cvalp, betap) = cvals_betas_5(&arrp);
        let a0 = D5[0] / (EPS5 + betap[0]).powi(2);
        let a1 = D5[1] / (EPS5 + betap[1]).powi(2);
        let a2 = D5[2] / (EPS5 + betap[2]).powi(2);
        let sp = a0 + a1 + a2;
        let ql = (a0 / sp) * cvalp[0] + (a1 / sp) * cvalp[1] + (a2 / sp) * cvalp[2];

        // Right reconstruction at i+1/2 (mirror) using WENO5-JS weights
        let arrm = [
            q[pid(i as isize + 3, n)],
            q[pid(i as isize + 2, n)],
            q[ip],
            q[i],
            q[pid(i as isize - 1, n)],
        ];
        let (cvalm, betam) = cvals_betas_5(&arrm);
        let a0m = D5[0] / (EPS5 + betam[0]).powi(2);
        let a1m = D5[1] / (EPS5 + betam[1]).powi(2);
        let a2m = D5[2] / (EPS5 + betam[2]).powi(2);
        let sm = a0m + a1m + a2m;
        // mirrored weights order
        let qr = (a2m / sm) * cvalm[0] + (a1m / sm) * cvalm[1] + (a0m / sm) * cvalm[2];

        let qu = if u >= 0.0 { ql } else { qr };
        fh[i] = u * qu;
    }
    fh
}

impl Scheme for Weno5Js {
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
            let fh = weno5_fd_flux_faces_1d_upwind(&qq, &uu);
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
            let gh = weno5_fd_flux_faces_1d_upwind(&qq, &vv);
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
