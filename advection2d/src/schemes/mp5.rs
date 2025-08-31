use crate::schemes::Scheme;
use crate::utils::{idx, pid};

/// Monotonicity preserving fifth-order scheme (MP5).
pub struct Mp5;

fn mp5_face_pos(q: &[f64], i: usize, n: usize) -> f64 {
    let im2 = q[pid(i as isize - 2, n)];
    let im1 = q[pid(i as isize - 1, n)];
    let i0 = q[i];
    let ip1 = q[pid(i as isize + 1, n)];
    let ip2 = q[pid(i as isize + 2, n)];
    let mut qf = (2.0 * im2 - 13.0 * im1 + 47.0 * i0 + 27.0 * ip1 - 3.0 * ip2) / 60.0;
    let min = im1.min(i0.min(ip1));
    let max = im1.max(i0.max(ip1));
    if qf < min {
        qf = min;
    } else if qf > max {
        qf = max;
    }
    qf
}

fn mp5_face_neg(q: &[f64], i: usize, n: usize) -> f64 {
    let ip2 = q[pid(i as isize + 2, n)];
    let ip1 = q[pid(i as isize + 1, n)];
    let i0 = q[i];
    let im1 = q[pid(i as isize - 1, n)];
    let im2 = q[pid(i as isize - 2, n)];
    let mut qf = (2.0 * ip2 - 13.0 * ip1 + 47.0 * i0 + 27.0 * im1 - 3.0 * im2) / 60.0;
    let min = im1.min(i0.min(ip1));
    let max = im1.max(i0.max(ip1));
    if qf < min {
        qf = min;
    } else if qf > max {
        qf = max;
    }
    qf
}

fn flux_x(q: &[f64], u: &[f64], _dx: f64, nx: usize, j: usize) -> Vec<f64> {
    let mut f = vec![0.0; nx];
    let mut qq = vec![0.0; nx];
    for i in 0..nx {
        qq[i] = q[idx(i, j, nx)];
    }
    for i in 0..nx {
        let ip = pid(i as isize + 1, nx);
        let idx0 = idx(i, j, nx);
        let idxp = idx(ip, j, nx);
        let u_face = (1.0 / 2.0) * (u[idx0] + u[idxp]);
        let qf = if u_face >= 0.0 {
            mp5_face_pos(&qq, i, nx)
        } else {
            mp5_face_neg(&qq, ip, nx)
        };
        f[i] = u_face * qf;
    }
    f
}

fn flux_y(q: &[f64], v: &[f64], _dy: f64, nx: usize, ny: usize, i: usize) -> Vec<f64> {
    let mut g = vec![0.0; ny];
    let mut qq = vec![0.0; ny];
    for j in 0..ny {
        qq[j] = q[idx(i, j, nx)];
    }
    for j in 0..ny {
        let jp = pid(j as isize + 1, ny);
        let idx0 = idx(i, j, nx);
        let idxp = idx(i, jp, nx);
        let v_face = (1.0 / 2.0) * (v[idx0] + v[idxp]);
        let qf = if v_face >= 0.0 {
            mp5_face_pos(&qq, j, ny)
        } else {
            mp5_face_neg(&qq, jp, ny)
        };
        g[j] = v_face * qf;
    }
    g
}

impl Scheme for Mp5 {
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
            let fx = flux_x(q, u, dx, nx, j);
            for i in 0..nx {
                let im = pid(i as isize - 1, nx);
                dfx[idx(i, j, nx)] = (fx[i] - fx[im]) / dx;
            }
        }
        let mut dfy = vec![0.0; nx * ny];
        for i in 0..nx {
            let fy = flux_y(q, v, dy, nx, ny, i);
            for j in 0..ny {
                let jm = pid(j as isize - 1, ny);
                dfy[idx(i, j, nx)] = (fy[j] - fy[jm]) / dy;
            }
        }
        for k in 0..nx * ny {
            out[k] = -(dfx[k] + dfy[k]);
        }
    }
}
