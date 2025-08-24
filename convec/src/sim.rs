use crate::config::{Config, SchemeType, VelocityCfg};
use crate::render::FrameWriter;
use crate::schemes::{Centered8, Scheme, Weno5Js};
use crate::shapes::init_field;
use crate::utils::idx;
use anyhow::Result;

#[derive(Debug, Clone, Copy)]
pub struct RunStats {
    pub qmin: f64,
    pub qmax: f64,
    pub mass0: f64,
    pub mass1: f64,
    pub l2: f64,
}

/// 2 次元スカラー移流方程式を解くメインループ。
/// 時間積分には Shu & Osher による 3 段の TVD Runge–Kutta 法\[1\]を用いる。
///
/// [1] C.-W. Shu and S. Osher, "Efficient implementation of essentially non-oscillatory
/// shock-capturing schemes, II", *Journal of Computational Physics*, 83(1), 32-78, 1989.
pub fn run(cfg: Config) -> Result<RunStats> {
    let nx = cfg.simulation.nx;
    let ny = cfg.simulation.ny;
    let lx = cfg.simulation.lx;
    let ly = cfg.simulation.ly;
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;

    let mut x = vec![0.0; nx];
    let mut y = vec![0.0; ny];
    for i in 0..nx {
        x[i] = (i as f64 + 0.5) * dx;
    }
    for j in 0..ny {
        y[j] = (j as f64 + 0.5) * dy;
    }

    let mut u = vec![0.0; nx * ny];
    let mut v = vec![0.0; nx * ny];
    match cfg.simulation.velocity {
        VelocityCfg::SolidRotation {
            omega,
            center_x,
            center_y,
        } => {
            for j in 0..ny {
                for i in 0..nx {
                    let k = idx(i, j, nx);
                    u[k] = -omega * (y[j] - center_y);
                    v[k] = omega * (x[i] - center_x);
                }
            }
        }
    }

    let mut q = init_field(&cfg.initial_condition, nx, ny, lx, ly);
    let q0 = q.clone();

    let mut umax = 0.0;
    for k in 0..nx * ny {
        let s = (u[k] * u[k] + v[k] * v[k]).sqrt();
        if s > umax {
            umax = s;
        }
    }
    let mut dt = cfg.simulation.cfl * dx.min(dy) / (umax + 1e-15);
    let t_end = cfg.simulation.rotations as f64;
    let steps = (t_end / dt).ceil() as usize;
    dt = t_end / steps as f64;
    println!("[info] steps={} dt={:.3e}", steps, dt);

    let scheme_box: Box<dyn Scheme> = match cfg.scheme.r#type {
        SchemeType::Centered8 => Box::new(Centered8),
        SchemeType::Weno5 => Box::new(Weno5Js),
    };

    let mut writer = FrameWriter::new(cfg.output.clone(), nx, ny)?;
    if writer.cfg.enable {
        writer.write_frame(&q, 0)?;
    }

    let mut rhs = vec![0.0; nx * ny];
    let mut q1 = vec![0.0; nx * ny];
    let mut q2 = vec![0.0; nx * ny];

    for n in 1..=steps {
        scheme_box.rhs(&q, &u, &v, dx, dy, nx, ny, &mut rhs);
        for k in 0..nx * ny {
            q1[k] = q[k] + dt * rhs[k];
        }

        scheme_box.rhs(&q1, &u, &v, dx, dy, nx, ny, &mut rhs);
        for k in 0..nx * ny {
            q2[k] = 0.75 * q[k] + 0.25 * (q1[k] + dt * rhs[k]);
        }

        scheme_box.rhs(&q2, &u, &v, dx, dy, nx, ny, &mut rhs);
        for k in 0..nx * ny {
            q[k] = (1.0 / 3.0) * q[k] + (2.0 / 3.0) * (q2[k] + dt * rhs[k]);
        }

        writer.maybe_write(&q, n)?;
    }

    let cell = dx * dy;
    let mut qmin = q[0];
    let mut qmax = q[0];
    let mut m0 = 0.0;
    let mut m1 = 0.0;
    let mut l2 = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j, nx);
            m0 += q0[k] * cell;
            m1 += q[k] * cell;
            if q[k] < qmin {
                qmin = q[k];
            }
            if q[k] > qmax {
                qmax = q[k];
            }
            let d = q[k] - q0[k];
            l2 += d * d * cell;
        }
    }
    let stats = RunStats {
        qmin,
        qmax,
        mass0: m0,
        mass1: m1,
        l2: l2.sqrt(),
    };
    println!("# results");
    println!("min={:.6e} max={:.6e}", stats.qmin, stats.qmax);
    println!("mass0={:.6e} mass1={:.6e} mass_err={:.6e}", stats.mass0, stats.mass1, stats.mass1 - stats.mass0);
    println!("L2 = {:.6e}", stats.l2);
    Ok(stats)
}
