//! 9次精度の加重基本的非振動 (WENO) スキーム。
//!
//! 実装方針:
//! - 空間再構成は WENO-Z（Borges ら）の重みを用いた 5 候補の 9点再構成。
//! - 数値流束は「面上の上流値」から構成する。面速度
//!   `u_{i+1/2} = 0.5(u_i + u_{i+1})` の符号で左右の再構成値（左: q_L, 右: q_R）
//!   を選び，`f_{i+1/2} = u_{i+1/2} * q_upwind` とする。
//!   これにより平滑場での不要な LLF 型分割の拡散を避ける。
//! - 周期境界は `pid` で処理し，インデックスの折り返しを明示する。
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

/// 1 次元列（または行）に対して，面上の上流値から数値流束を組み立てる。
fn weno9z_fd_flux_faces_1d_upwind(q: &[f64], vel: &[f64]) -> Vec<f64> {
    let n = q.len();
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let ip = pid(i as isize + 1, n);
        let u = 0.5 * (vel[i] + vel[ip]);
        let mut arrp = [0.0; 9];
        for k in 0..9 { arrp[k] = q[pid(i as isize + k as isize - 4, n)]; }
        let ql = reconstruct(&arrp);
        let mut arrm = [0.0; 9];
        for k in 0..9 { arrm[k] = q[pid(i as isize + 5 - k as isize, n)]; }
        let qr = reconstruct(&arrm);
        let qu = if u >= 0.0 { ql } else { qr };
        fh[i] = u * qu;
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
            // x 方向の 1D 列に投影して上流面流束を形成
            let mut qq = vec![0.0; nx];
            let mut uu = vec![0.0; nx];
            for i in 0..nx {
                let k = idx(i, j, nx);
                qq[i] = q[k];
                uu[i] = u[k];
            }
            let fh = weno9z_fd_flux_faces_1d_upwind(&qq, &uu);
            for i in 0..nx {
                let im = pid(i as isize - 1, nx);
                dfx[idx(i, j, nx)] = (fh[i] - fh[im]) / dx;
            }
        }
        let mut dfy = vec![0.0; nx * ny];
        for i in 0..nx {
            // y 方向の 1D 列に投影して上流面流束を形成
            let mut qq = vec![0.0; ny];
            let mut vv = vec![0.0; ny];
            for j in 0..ny {
                let k = idx(i, j, nx);
                qq[j] = q[k];
                vv[j] = v[k];
            }
            let gh = weno9z_fd_flux_faces_1d_upwind(&qq, &vv);
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
