use crate::schemes::Scheme;
use crate::utils::{idx, pid};

/// Compute cubic interpolated face value using cell-centered values and derivatives.
fn cip_face(ql: f64, qr: f64, dql: f64, dqr: f64, dx: f64) -> f64 {
    let a = (3.0 * (qr - ql) / dx - (2.0 * dql + dqr)) / dx;
    let b = (2.0 * (ql - qr) / dx + (dql + dqr)) / (dx * dx);
    let x = 0.5 * dx;
    ql + dql * x + a * x * x + b * x * x * x
}

/// Limit left and right slopes so that the cubic Hermite profile remains monotonic.
fn monotonic_slopes(ql: f64, qr: f64, mut dl: f64, mut dr: f64, dx: f64) -> (f64, f64) {
    let delta = (qr - ql) / dx;
    if delta == 0.0 {
        return (0.0, 0.0);
    }
    let mut r = dl / delta;
    let mut s = dr / delta;
    if r * s <= 0.0 {
        return (0.0, 0.0);
    }
    let sum = r + s;
    if sum > 3.0 {
        let k = 3.0 / sum;
        r *= k;
        s *= k;
    }
    dl = r * delta;
    dr = s * delta;
    (dl, dr)
}

fn flux_x(
    q: &[f64],
    dqx: &[f64],
    u: &[f64],
    dx: f64,
    nx: usize,
    j: usize,
    limiter: bool,
) -> Vec<f64> {
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

fn flux_y(
    q: &[f64],
    dqy: &[f64],
    v: &[f64],
    dy: f64,
    nx: usize,
    ny: usize,
    i: usize,
    limiter: bool,
) -> Vec<f64> {
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

/// Sixth-order central difference of the x-derivative.
fn derivative_x6(q: &[f64], dx: f64, nx: usize, ny: usize) -> Vec<f64> {
    let mut dqx = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let im3 = pid(i as isize - 3, nx);
            let im2 = pid(i as isize - 2, nx);
            let im1 = pid(i as isize - 1, nx);
            let ip1 = pid(i as isize + 1, nx);
            let ip2 = pid(i as isize + 2, nx);
            let ip3 = pid(i as isize + 3, nx);
            let k = idx(i, j, nx);
            dqx[k] = (
                q[idx(im3, j, nx)]
                    - 9.0 * q[idx(im2, j, nx)]
                    + 45.0 * q[idx(im1, j, nx)]
                    - 45.0 * q[idx(ip1, j, nx)]
                    + 9.0 * q[idx(ip2, j, nx)]
                    - q[idx(ip3, j, nx)]
            ) / (60.0 * dx);
        }
    }
    dqx
}

/// Sixth-order central difference of the y-derivative.
fn derivative_y6(q: &[f64], dy: f64, nx: usize, ny: usize) -> Vec<f64> {
    let mut dqy = vec![0.0; nx * ny];
    for i in 0..nx {
        for j in 0..ny {
            let jm3 = pid(j as isize - 3, ny);
            let jm2 = pid(j as isize - 2, ny);
            let jm1 = pid(j as isize - 1, ny);
            let jp1 = pid(j as isize + 1, ny);
            let jp2 = pid(j as isize + 2, ny);
            let jp3 = pid(j as isize + 3, ny);
            let k = idx(i, j, nx);
            dqy[k] = (
                q[idx(i, jm3, nx)]
                    - 9.0 * q[idx(i, jm2, nx)]
                    + 45.0 * q[idx(i, jm1, nx)]
                    - 45.0 * q[idx(i, jp1, nx)]
                    + 9.0 * q[idx(i, jp2, nx)]
                    - q[idx(i, jp3, nx)]
            ) / (60.0 * dy);
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
/// Higher-order conservative CIP-CSL2 scheme.
pub struct CipCsl2;
/// CIP-CSL2 scheme with monotonicity-constrained Hermite interpolation.
pub struct CipCsl2Mh;

fn cip_rhs(
    q: &[f64],
    u: &[f64],
    v: &[f64],
    dx: f64,
    dy: f64,
    nx: usize,
    ny: usize,
    limiter: bool,
    out: &mut [f64],
) {
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

fn cip_csl2_rhs(
    q: &[f64],
    u: &[f64],
    v: &[f64],
    dx: f64,
    dy: f64,
    nx: usize,
    ny: usize,
    monotonic: bool,
    out: &mut [f64],
) {
    let dqx = derivative_x6(q, dx, nx, ny);
    let dqy = derivative_y6(q, dy, nx, ny);
    let mut dfx = vec![0.0; nx * ny];
    for j in 0..ny {
        let mut f = vec![0.0; nx];
        for i in 0..nx {
            let ip = pid(i as isize + 1, nx);
            let idx0 = idx(i, j, nx);
            let idxp = idx(ip, j, nx);
            let u_face = 0.5 * (u[idx0] + u[idxp]);
            let (ql, qr, mut dql, mut dqr) = if u_face >= 0.0 {
                (q[idx0], q[idxp], dqx[idx0], dqx[idxp])
            } else {
                (q[idxp], q[idx0], dqx[idxp], dqx[idx0])
            };
            if monotonic {
                let (dl, dr) = monotonic_slopes(ql, qr, dql, dqr, dx);
                dql = dl;
                dqr = dr;
            }
            let mut qf = cip_face(ql, qr, dql, dqr, dx);
            if monotonic {
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
        for i in 0..nx {
            let im = pid(i as isize - 1, nx);
            dfx[idx(i, j, nx)] = (f[i] - f[im]) / dx;
        }
    }
    let mut dfy = vec![0.0; nx * ny];
    for i in 0..nx {
        let mut g = vec![0.0; ny];
        for j in 0..ny {
            let jp = pid(j as isize + 1, ny);
            let idx0 = idx(i, j, nx);
            let idxp = idx(i, jp, nx);
            let v_face = 0.5 * (v[idx0] + v[idxp]);
            let (ql, qr, mut dql, mut dqr) = if v_face >= 0.0 {
                (q[idx0], q[idxp], dqy[idx0], dqy[idxp])
            } else {
                (q[idxp], q[idx0], dqy[idxp], dqy[idx0])
            };
            if monotonic {
                let (dl, dr) = monotonic_slopes(ql, qr, dql, dqr, dy);
                dql = dl;
                dqr = dr;
            }
            let mut qf = cip_face(ql, qr, dql, dqr, dy);
            if monotonic {
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
        for j in 0..ny {
            let jm = pid(j as isize - 1, ny);
            dfy[idx(i, j, nx)] = (g[j] - g[jm]) / dy;
        }
    }
    for k in 0..nx * ny {
        out[k] = -(dfx[k] + dfy[k]);
    }
}

impl Scheme for Cip {
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
        cip_rhs(q, u, v, dx, dy, nx, ny, false, out);
    }
}

impl Scheme for CipCsl {
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
        cip_rhs(q, u, v, dx, dy, nx, ny, false, out);
    }
}

impl Scheme for CipB {
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
        cip_rhs(q, u, v, dx, dy, nx, ny, true, out);
    }
}

impl Scheme for CipCsl2 {
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
        cip_csl2_rhs(q, u, v, dx, dy, nx, ny, false, out);
    }
}

impl Scheme for CipCsl2Mh {
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
        cip_csl2_rhs(q, u, v, dx, dy, nx, ny, true, out);
    }
}
