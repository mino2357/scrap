use crate::schemes::Scheme;
use crate::utils::{idx, pid};

/// Compute cubic interpolated face value using cell-centered values and derivatives.
fn cip_face(ql: f64, qr: f64, dql: f64, dqr: f64, dx: f64) -> f64 {
    let a = (3.0 * (qr - ql) / dx - (2.0 * dql + dqr)) / dx;
    let b = (2.0 * (ql - qr) / dx + (dql + dqr)) / (dx * dx);
    let x = 0.5 * dx;
    ql + dql * x + a * x * x + b * x * x * x
}

fn flux_x(q: &[f64], dqx: &[f64], u: &[f64], dx: f64, nx: usize, j: usize, limiter: bool) -> Vec<f64> {
    let mut f = vec![0.0; nx];
    for i in 0..nx {
        let ip = pid(i as isize + 1, nx);
        let idx0 = idx(i, j, nx);
        let idxp = idx(ip, j, nx);
        let u_face = 0.5 * (u[idx0] + u[idxp]);
        let (ql, qr, dql, dqr) = if u_face >= 0.0 {
            (q[idx0], q[idxp], dqx[idx0], dqx[idxp])
        } else {
            (q[idxp], q[idx0], dqx[idxp], dqx[idx0])
        };
        let mut qf = cip_face(ql, qr, dql, dqr, dx);
        if limiter {
            let min = ql.min(qr);
            let max = ql.max(qr);
            if qf < min {
                qf = min;
            } else if qf > max {
                qf = max;
            }
        }
        f[i] = u_face * qf;
    }
    f
}

fn flux_y(q: &[f64], dqy: &[f64], v: &[f64], dy: f64, nx: usize, ny: usize, i: usize, limiter: bool) -> Vec<f64> {
    let mut g = vec![0.0; ny];
    for j in 0..ny {
        let jp = pid(j as isize + 1, ny);
        let idx0 = idx(i, j, nx);
        let idxp = idx(i, jp, nx);
        let v_face = 0.5 * (v[idx0] + v[idxp]);
        let (ql, qr, dql, dqr) = if v_face >= 0.0 {
            (q[idx0], q[idxp], dqy[idx0], dqy[idxp])
        } else {
            (q[idxp], q[idx0], dqy[idxp], dqy[idx0])
        };
        let mut qf = cip_face(ql, qr, dql, dqr, dy);
        if limiter {
            let min = ql.min(qr);
            let max = ql.max(qr);
            if qf < min {
                qf = min;
            } else if qf > max {
                qf = max;
            }
        }
        g[j] = v_face * qf;
    }
    g
}

fn derivative_x(q: &[f64], dx: f64, nx: usize, ny: usize) -> Vec<f64> {
    let mut dqx = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let im = pid(i as isize - 1, nx);
            let ip = pid(i as isize + 1, nx);
            let k = idx(i, j, nx);
            dqx[k] = (q[idx(ip, j, nx)] - q[idx(im, j, nx)]) / (2.0 * dx);
        }
    }
    dqx
}

fn derivative_y(q: &[f64], dy: f64, nx: usize, ny: usize) -> Vec<f64> {
    let mut dqy = vec![0.0; nx * ny];
    for j in 0..ny {
        let jm = pid(j as isize - 1, ny);
        let jp = pid(j as isize + 1, ny);
        for i in 0..nx {
            let k = idx(i, j, nx);
            dqy[k] = (q[idx(i, jp, nx)] - q[idx(i, jm, nx)]) / (2.0 * dy);
        }
    }
    dqy
}

/// Cubic interpolated propagation scheme.
pub struct Cip;
/// Conservative CIP scheme preserving cell averages.
pub struct CipCsl;
/// Bounded CIP scheme with a min–max limiter.
pub struct CipB;

fn cip_rhs(q: &[f64], u: &[f64], v: &[f64], dx: f64, dy: f64, nx: usize, ny: usize, limiter: bool, out: &mut [f64]) {
    let dqx = derivative_x(q, dx, nx, ny);
    let dqy = derivative_y(q, dy, nx, ny);
    let mut dfx = vec![0.0; nx * ny];
    for j in 0..ny {
        let fx = flux_x(q, &dqx, u, dx, nx, j, limiter);
        for i in 0..nx {
            let im = pid(i as isize - 1, nx);
            dfx[idx(i, j, nx)] = (fx[i] - fx[im]) / dx;
        }
    }
    let mut dfy = vec![0.0; nx * ny];
    for i in 0..nx {
        let fy = flux_y(q, &dqy, v, dy, nx, ny, i, limiter);
        for j in 0..ny {
            let jm = pid(j as isize - 1, ny);
            dfy[idx(i, j, nx)] = (fy[j] - fy[jm]) / dy;
        }
    }
    for k in 0..nx * ny {
        out[k] = -(dfx[k] + dfy[k]);
    }
}

impl Scheme for Cip {
    fn rhs(&self, q: &[f64], u: &[f64], v: &[f64], dx: f64, dy: f64, nx: usize, ny: usize, out: &mut [f64]) {
        cip_rhs(q, u, v, dx, dy, nx, ny, false, out);
    }
}

impl Scheme for CipCsl {
    fn rhs(&self, q: &[f64], u: &[f64], v: &[f64], dx: f64, dy: f64, nx: usize, ny: usize, out: &mut [f64]) {
        cip_rhs(q, u, v, dx, dy, nx, ny, false, out);
    }
}

impl Scheme for CipB {
    fn rhs(&self, q: &[f64], u: &[f64], v: &[f64], dx: f64, dy: f64, nx: usize, ny: usize, out: &mut [f64]) {
        cip_rhs(q, u, v, dx, dy, nx, ny, true, out);
    }
}

