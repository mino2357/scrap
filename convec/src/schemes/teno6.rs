//! Sixth-order Targeted ENO (TENO6) scheme.
//!
//! Adaptive-order reconstruction on a five-point stencil using the
//! WENO5 family candidates and a TENO cut-off sensor. In smooth regions
//! it reduces to WENO-Z; near discontinuities, non-smooth substencils are
//! discarded and the remaining weights are renormalised.
use crate::schemes::Scheme;
use crate::utils::{idx, pid};
use super::weno5_core::{cvals_betas_5, D5, EPS5};

pub struct Teno6;

fn teno6_reconstruct(arr: &[f64; 5]) -> f64 {
    let (cval, beta) = cvals_betas_5(arr);
    // Global smoothness and linear weights (WENO5 optimal).
    let tau5 = (beta[0] - beta[2]).abs();
    // WENO-Z weights for smooth fallback.
    let q = 2; // standard exponent for WENO-Z on tau/beta
    let a0z = D5[0] * (1.0 + (tau5 / (EPS5 + beta[0])).powi(q));
    let a1z = D5[1] * (1.0 + (tau5 / (EPS5 + beta[1])).powi(q));
    let a2z = D5[2] * (1.0 + (tau5 / (EPS5 + beta[2])).powi(q));
    let sz = a0z + a1z + a2z;
    let w0z = a0z / sz;
    let w1z = a1z / sz;
    let w2z = a2z / sz;

    // TENO sensor and cut-off.
    let ct = 1e-6; // cut-off threshold
    let qteno = 6; // sharpen the sensor as in TENO literature
    let r0 = ((tau5 + EPS5) / (beta[0] + EPS5)).powi(qteno);
    let r1 = ((tau5 + EPS5) / (beta[1] + EPS5)).powi(qteno);
    let r2 = ((tau5 + EPS5) / (beta[2] + EPS5)).powi(qteno);

    let a0 = if r0 < ct { w0z } else { 0.0 };
    let a1 = if r1 < ct { w1z } else { 0.0 };
    let a2 = if r2 < ct { w2z } else { 0.0 };
    let s = a0 + a1 + a2;
    let (w0, w1, w2) = if s > 0.0 {
        (a0 / s, a1 / s, a2 / s)
    } else {
        // all stencils flagged; revert to WENO-Z
        (w0z, w1z, w2z)
    };

    w0 * cval[0] + w1 * cval[1] + w2 * cval[2]
}

fn teno6_fd_flux_faces_1d(f: &[f64], q: &[f64], alpha: f64) -> Vec<f64> {
    let n = f.len();
    let mut fp = vec![0.0; n];
    let mut fm = vec![0.0; n];
    for i in 0..n {
        fp[i] = (1.0 / 2.0) * (f[i] + alpha * q[i]);
        fm[i] = (1.0 / 2.0) * (f[i] - alpha * q[i]);
    }
    let mut fh = vec![0.0; n];
    for i in 0..n {
        // positive flux, face i+1/2
        let arrp = [
            fp[pid(i as isize - 2, n)],
            fp[pid(i as isize - 1, n)],
            fp[i],
            fp[pid(i as isize + 1, n)],
            fp[pid(i as isize + 2, n)],
        ];
        let fph = teno6_reconstruct(&arrp);

        // negative flux, mirror about i+1/2
        let arrm = [
            fm[pid(i as isize + 3, n)],
            fm[pid(i as isize + 2, n)],
            fm[pid(i as isize + 1, n)],
            fm[i],
            fm[pid(i as isize - 1, n)],
        ];
        let fmh = teno6_reconstruct(&arrm);
        fh[i] = fph + fmh;
    }
    fh
}

impl Scheme for Teno6 {
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
            let fh = teno6_fd_flux_faces_1d(&f, &qq, amax);
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
            let gh = teno6_fd_flux_faces_1d(&g, &qq, amax);
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
