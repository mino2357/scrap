//! 7次精度の加重基本的非振動 (WENO) スキームを用いた移流項評価。
//! Borges らによる WENO-Z 重み\[1\]を用いて滑らかな領域での精度を改善する。
//!
//! [1] R. Borges, M. Carmona, B. Costa, and W. S. Don,
//!     "An improved weighted essentially non-oscillatory scheme for
//!     hyperbolic conservation laws", *Journal of Computational Physics*,
//!     227(6), 3191-3211, 2008.
use crate::schemes::Scheme;
use crate::utils::{idx, pid};

pub struct Weno7Z;

fn beta(a1: f64, a2: f64, a3: f64) -> f64 {
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
    let c0 = (-5.0 * im3 + 21.0 * im2 - 35.0 * im1 + 35.0 * i0) / 16.0;
    let c1 = (im2 - 5.0 * im1 + 15.0 * i0 + 5.0 * ip1) / 16.0;
    let c2 = (-im1 + 9.0 * i0 + 9.0 * ip1 - ip2) / 16.0;
    let c3 = (5.0 * i0 + 15.0 * ip1 - 5.0 * ip2 + ip3) / 16.0;
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
    let eps = 1e-6;
    let tau7 = (b0 - b3).abs();
    let d0 = 1.0 / 35.0;
    let d1 = 12.0 / 35.0;
    let d2 = 18.0 / 35.0;
    let d3 = 4.0 / 35.0;
    let a0 = d0 * (1.0 + (tau7 / (b0 + eps)).powi(2));
    let a1 = d1 * (1.0 + (tau7 / (b1 + eps)).powi(2));
    let a2 = d2 * (1.0 + (tau7 / (b2 + eps)).powi(2));
    let a3 = d3 * (1.0 + (tau7 / (b3 + eps)).powi(2));
    let s = a0 + a1 + a2 + a3;
    let w0 = a0 / s;
    let w1 = a1 / s;
    let w2 = a2 / s;
    let w3 = a3 / s;
    w0 * c0 + w1 * c1 + w2 * c2 + w3 * c3
}

fn weno7z_fd_flux_faces_1d(f: &[f64], q: &[f64], alpha: f64) -> Vec<f64> {
    let n = f.len();
    let mut fp = vec![0.0; n];
    let mut fm = vec![0.0; n];
    for i in 0..n {
        fp[i] = (1.0 / 2.0) * (f[i] + alpha * q[i]);
        fm[i] = (1.0 / 2.0) * (f[i] - alpha * q[i]);
    }
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let arrp = [
            fp[pid(i as isize - 3, n)],
            fp[pid(i as isize - 2, n)],
            fp[pid(i as isize - 1, n)],
            fp[i],
            fp[pid(i as isize + 1, n)],
            fp[pid(i as isize + 2, n)],
            fp[pid(i as isize + 3, n)],
        ];
        let fph = reconstruct(&arrp);
        let arrm = [
            fm[pid(i as isize + 4, n)],
            fm[pid(i as isize + 3, n)],
            fm[pid(i as isize + 2, n)],
            fm[pid(i as isize + 1, n)],
            fm[i],
            fm[pid(i as isize - 1, n)],
            fm[pid(i as isize - 2, n)],
        ];
        let fmh = reconstruct(&arrm);
        fh[i] = fph + fmh;
    }
    fh
}

impl Scheme for Weno7Z {
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
            let fh = weno7z_fd_flux_faces_1d(&f, &qq, amax);
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
            let gh = weno7z_fd_flux_faces_1d(&g, &qq, amax);
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
