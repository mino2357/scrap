//! 5次精度の加重基本的非振動 (WENO) スキームを用いた移流項評価。
//! 実装は Jiang & Shu による古典的な WENO-JS スキーム\[1\]に基づいている。
//!
//! [1] G.-S. Jiang and C.-W. Shu, "Efficient implementation of Weighted ENO schemes",
//! *Journal of Computational Physics*, 126(1), 202-228, 1996.
use crate::schemes::Scheme;
use crate::utils::{idx, pid};

pub struct Weno5Js;

fn weno5_fd_flux_faces_1d(f: &[f64], q: &[f64], alpha: f64) -> Vec<f64> {
    let n = f.len();
    let eps = 1e-6;
    let mut fp = vec![0.0; n];
    let mut fm = vec![0.0; n];
    for i in 0..n {
        fp[i] = 0.5 * (f[i] + alpha * q[i]);
        fm[i] = 0.5 * (f[i] - alpha * q[i]);
    }
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let im2 = fp[pid(i as isize - 2, n)];
        let im1 = fp[pid(i as isize - 1, n)];
        let i0 = fp[i];
        let ip1 = fp[pid(i as isize + 1, n)];
        let ip2 = fp[pid(i as isize + 2, n)];
        let c0 = (2.0 * im2 - 7.0 * im1 + 11.0 * i0) / 6.0;
        let c1 = (-im1 + 5.0 * i0 + 2.0 * ip1) / 6.0;
        let c2 = (2.0 * i0 + 5.0 * ip1 - ip2) / 6.0;
        let b0 = (13.0 / 12.0) * (im2 - 2.0 * im1 + i0).powi(2)
            + 0.25 * (im2 - 4.0 * im1 + 3.0 * i0).powi(2);
        let b1 = (13.0 / 12.0) * (im1 - 2.0 * i0 + ip1).powi(2) + 0.25 * (im1 - ip1).powi(2);
        let b2 = (13.0 / 12.0) * (i0 - 2.0 * ip1 + ip2).powi(2)
            + 0.25 * (3.0 * i0 - 4.0 * ip1 + ip2).powi(2);
        let a0 = 0.1 / (eps + b0).powi(2);
        let a1 = 0.6 / (eps + b1).powi(2);
        let a2 = 0.3 / (eps + b2).powi(2);
        let s = a0 + a1 + a2;
        let w0 = a0 / s;
        let w1 = a1 / s;
        let w2 = a2 / s;
        let fph = w0 * c0 + w1 * c1 + w2 * c2;

        let im1m = fm[pid(i as isize - 1, n)];
        let i0m = fm[i];
        let ip1m = fm[pid(i as isize + 1, n)];
        let ip2m = fm[pid(i as isize + 2, n)];
        let ip3m = fm[pid(i as isize + 3, n)];
        let b0m = (13.0 / 12.0) * (i0m - 2.0 * ip1m + ip2m).powi(2)
            + 0.25 * (3.0 * i0m - 4.0 * ip1m + ip2m).powi(2);
        let b1m = (13.0 / 12.0) * (im1m - 2.0 * i0m + ip1m).powi(2) + 0.25 * (im1m - ip1m).powi(2);
        let b2m = (13.0 / 12.0) * (im1m - 2.0 * ip1m + ip2m).powi(2)
            + 0.25 * (im1m - 4.0 * i0m + 3.0 * ip1m).powi(2);
        let a0m = 0.1 / (eps + b0m).powi(2);
        let a1m = 0.6 / (eps + b1m).powi(2);
        let a2m = 0.3 / (eps + b2m).powi(2);
        let sm = a0m + a1m + a2m;
        let w0m = a0m / sm;
        let w1m = a1m / sm;
        let w2m = a2m / sm;
        let cm0 = (-im1m + 5.0 * i0m + 2.0 * ip1m) / 6.0;
        let cm1 = (2.0 * i0m + 5.0 * ip1m - ip2m) / 6.0;
        let cm2 = (11.0 * ip1m - 7.0 * ip2m + 2.0 * ip3m) / 6.0;
        let fmh = w2m * cm0 + w1m * cm1 + w0m * cm2;
        fh[i] = fph + fmh;
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
            let fh = weno5_fd_flux_faces_1d(&f, &qq, amax);
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
            let gh = weno5_fd_flux_faces_1d(&g, &qq, amax);
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
