use crate::schemes::Scheme;
use crate::utils::{idx, pid};

#[inline]
fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 { 0.0 } else { a.signum() * a.abs().min(b.abs()) }
}

#[inline]
fn minmod3(a: f64, b: f64, c: f64) -> f64 {
    minmod(minmod(a, b), c)
}

#[inline]
fn slope_minmod(q: &[f64], i: usize, n: usize) -> f64 {
    let im1 = pid(i as isize - 1, n);
    let ip1 = pid(i as isize + 1, n);
    let d_m = q[i] - q[im1];
    let d_p = q[ip1] - q[i];
    let d_c = 0.5 * (q[ip1] - q[im1]);
    minmod3(d_m, d_c, d_p)
}

#[inline]
fn slope_van_leer(q: &[f64], i: usize, n: usize) -> f64 {
    // van Leer harmonic mean limiter (symmetric in +/- directions)
    let im1 = pid(i as isize - 1, n);
    let ip1 = pid(i as isize + 1, n);
    let a = q[i] - q[im1];
    let b = q[ip1] - q[i];
    if a * b <= 0.0 { 0.0 } else { (2.0 * a * b) / (a + b) }
}

fn tvd_flux_faces_1d(q: &[f64], u: &[f64], use_vanleer: bool) -> Vec<f64> {
    let n = q.len();
    let mut fh = vec![0.0; n];
    // Precompute node slopes (symmetric reconstruction)
    let mut s = vec![0.0; n];
    for i in 0..n {
        s[i] = if use_vanleer { slope_van_leer(q, i, n) } else { slope_minmod(q, i, n) };
    }
    for i in 0..n {
        let ip1 = pid(i as isize + 1, n);
        let u_face = 0.5 * (u[i] + u[ip1]);
        let ql = q[i] + 0.5 * s[i];
        let qr = q[ip1] - 0.5 * s[ip1];
        fh[i] = if u_face >= 0.0 { u_face * ql } else { u_face * qr };
    }
    fh
}

fn tvd_rhs(
    q: &[f64],
    u: &[f64],
    v: &[f64],
    dx: f64,
    dy: f64,
    nx: usize,
    ny: usize,
    out: &mut [f64],
    use_vanleer: bool,
) {
        let mut dfx = vec![0.0; nx * ny];
        for j in 0..ny {
            let mut qrow = vec![0.0; nx];
            let mut urow = vec![0.0; nx];
            for i in 0..nx {
                let k = idx(i, j, nx);
                qrow[i] = q[k];
                urow[i] = u[k];
            }
        let fh = tvd_flux_faces_1d(&qrow, &urow, use_vanleer);
            for i in 0..nx {
                let im = pid(i as isize - 1, nx);
                dfx[idx(i, j, nx)] = (fh[i] - fh[im]) / dx;
            }
        }
    let mut dfy = vec![0.0; nx * ny];
        for i in 0..nx {
            let mut qcol = vec![0.0; ny];
            let mut vcol = vec![0.0; ny];
            for j in 0..ny {
                let k = idx(i, j, nx);
                qcol[j] = q[k];
                vcol[j] = v[k];
            }
            let gh = tvd_flux_faces_1d(&qcol, &vcol, use_vanleer);
            for j in 0..ny {
                let jm = pid(j as isize - 1, ny);
                dfy[idx(i, j, nx)] = (gh[j] - gh[jm]) / dy;
            }
        }
    for k in 0..nx * ny {
        out[k] = -(dfx[k] + dfy[k]);
    }
}

/// minmod リミッターを用いたTVDスキーム。
pub struct TvdMinmod;

impl Scheme for TvdMinmod {
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
        tvd_rhs(q, u, v, dx, dy, nx, ny, out, false);
    }
}

/// van Leer リミッターを用いたTVDスキーム。
pub struct TvdVanLeer;

impl Scheme for TvdVanLeer {
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
        tvd_rhs(q, u, v, dx, dy, nx, ny, out, true);
    }
}
