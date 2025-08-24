use crate::schemes::Scheme;
use crate::utils::{idx, pid};

fn limiter_minmod(r: f64) -> f64 {
    r.max(0.0).min(1.0)
}

fn limiter_van_leer(r: f64) -> f64 {
    if r <= 0.0 {
        0.0
    } else {
        (r + r.abs()) / (1.0 + r.abs())
    }
}

fn tvd_flux_faces_1d(q: &[f64], u: &[f64], limiter: fn(f64) -> f64) -> Vec<f64> {
    let n = q.len();
    let eps = 1e-12;
    let mut fh = vec![0.0; n];
    for i in 0..n {
        let im1 = pid(i as isize - 1, n);
        let ip1 = pid(i as isize + 1, n);
        let ip2 = pid(i as isize + 2, n);
        let dqf = q[ip1] - q[i];
        let dqb = q[i] - q[im1];
        let r = if dqf.abs() < eps { 0.0 } else { dqb / dqf };
        let ql = q[i] + 0.5 * limiter(r) * dqf;
        let dqf1 = q[ip2] - q[ip1];
        let dqb1 = q[ip1] - q[i];
        let r1 = if dqf1.abs() < eps { 0.0 } else { dqb1 / dqf1 };
        let qr = q[ip1] - 0.5 * limiter(r1) * (q[ip1] - q[i]);
        let u_face = 0.5 * (u[i] + u[ip1]);
        fh[i] = if u_face >= 0.0 {
            u_face * ql
        } else {
            u_face * qr
        };
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
    limiter: fn(f64) -> f64,
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
        let fh = tvd_flux_faces_1d(&qrow, &urow, limiter);
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
        let gh = tvd_flux_faces_1d(&qcol, &vcol, limiter);
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
        tvd_rhs(q, u, v, dx, dy, nx, ny, out, limiter_minmod);
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
        tvd_rhs(q, u, v, dx, dy, nx, ny, out, limiter_van_leer);
    }
}
