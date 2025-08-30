//! 9次精度の加重基本的非振動 (WENO) スキームを用いた移流項評価。
//! Borges らによる WENO-Z 重みを用いて滑らかな領域での精度を改善する。
use crate::schemes::Scheme;
use crate::utils::{idx, pid};
use super::weno9_core::{cvals_betas, D, EPS};

pub struct Weno9Z;

fn reconstruct(arr: &[f64; 9]) -> f64 {
    let (cval, beta) = cvals_betas(arr);
    let tau = (beta[0] - beta[4]).abs();
    let mut alpha = [0.0; 5];
    let mut sum = 0.0;
    for k in 0..5 {
        alpha[k] = D[k] * (1.0 + (tau / (beta[k] + EPS)).powi(2));
        sum += alpha[k];
    }
    let mut w = [0.0; 5];
    for k in 0..5 {
        w[k] = alpha[k] / sum;
    }
    w.iter().zip(cval.iter()).map(|(w, c)| w * c).sum()
}

fn weno9z_fd_flux_faces_1d(f: &[f64], q: &[f64], alpha: f64) -> Vec<f64> {
    let n = f.len();
    let mut fp = vec![0.0; n];
    let mut fm = vec![0.0; n];
    for i in 0..n {
        fp[i] = (1.0 / 2.0) * (f[i] + alpha * q[i]);
        fm[i] = (1.0 / 2.0) * (f[i] - alpha * q[i]);
    }
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let mut arrp = [0.0; 9];
        for k in 0..9 {
            arrp[k] = fp[pid(i as isize + k as isize - 4, n)];
        }
        let fph = reconstruct(&arrp);
        let mut arrm = [0.0; 9];
        for k in 0..9 {
            arrm[k] = fm[pid(i as isize + 5 - k as isize, n)];
        }
        let fmh = reconstruct(&arrm);
        fh[i] = fph + fmh;
    }
    fh
}

impl Scheme for Weno9Z {
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
            let fh = weno9z_fd_flux_faces_1d(&f, &qq, amax);
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
            let gh = weno9z_fd_flux_faces_1d(&g, &qq, amax);
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
