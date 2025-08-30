//! シミュレーションのメインルーチン。

use crate::config::{Config, SchemeType, TimeIntegrator, VelocityCfg};
use crate::render::FrameWriter;
use crate::schemes::{
    Centered6, Centered8, Centered10, Centered12, Centered14, Cip, CipB, CipCsl, CipCsl2,
    CipCsl2Mh, Mp5, Scheme, Teno6, Teno7A, Teno8A, Teno9A, TvdMinmod, TvdVanLeer, Upwind1,
    Upwind3x3, Weno5Js, Weno5Z, Weno7Z, Weno9Z,
};
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
/// 時間積分には強安定性保存 (SSP) Runge–Kutta 法を用いる。
/// 現在は 3 段 3 次の SSPRK(3,3) 法\[1\]のみをサポートする。
///
/// [1] C.-W. Shu and S. Osher, "Efficient implementation of essentially non-oscillatory
/// shock-capturing schemes, II", *Journal of Computational Physics*, 83(1), 32-78, 1989.
pub fn run(cfg: Config) -> Result<RunStats> {
    let nx = cfg.simulation.nx;
    let ny = cfg.simulation.ny;
    let lx = cfg.simulation.lx;
    let ly = cfg.simulation.ly;
    // 格子幅。正方形格子を仮定
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;

    // 格子中心座標を計算しておく。
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
    // 設定に従い速度場を生成。現在は剛体回転のみをサポート。
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

    // 速度の最大値から時間刻み幅を決定
    let mut umax = 0.0;
    for k in 0..nx * ny {
        let s = (u[k] * u[k] + v[k] * v[k]).sqrt();
        if s > umax {
            umax = s;
        }
    }
    let mut dt = cfg.simulation.cfl * dx.min(dy) / (umax + 1e-15);
    // 1 回転を単位時間として終了時刻を決定
    let t_end = cfg.simulation.rotations as f64;
    let steps = (t_end / dt).ceil() as usize;
    dt = t_end / steps as f64;
    println!("[info] steps={} dt={:.3e}", steps, dt);

    let scheme_box: Box<dyn Scheme> = match cfg.scheme.r#type {
        SchemeType::Centered6 => Box::new(Centered6),
        SchemeType::Centered8 => Box::new(Centered8),
        SchemeType::Centered10 => Box::new(Centered10),
        SchemeType::Centered12 => Box::new(Centered12),
        SchemeType::Centered14 => Box::new(Centered14),
        SchemeType::Weno5 => Box::new(Weno5Js),
        SchemeType::Weno5Z => Box::new(Weno5Z),
        SchemeType::Weno7Z => Box::new(Weno7Z),
        SchemeType::Weno9Z => Box::new(Weno9Z),
        SchemeType::Upwind1 => Box::new(Upwind1),
        SchemeType::Upwind3x3 => Box::new(Upwind3x3),
        SchemeType::Cip => Box::new(Cip),
        SchemeType::CipCsl => Box::new(CipCsl),
        SchemeType::CipB => Box::new(CipB),
        SchemeType::CipCsl2 => Box::new(CipCsl2),
        SchemeType::CipCsl2Mh => Box::new(CipCsl2Mh),
        SchemeType::Mp5 => Box::new(Mp5),
        SchemeType::TvdMinmod => Box::new(TvdMinmod),
        SchemeType::TvdVanLeer => Box::new(TvdVanLeer),
        SchemeType::Teno6 => Box::new(Teno6),
        SchemeType::Teno7A => Box::new(Teno7A),
        SchemeType::Teno8A => Box::new(Teno8A),
        SchemeType::Teno9A => Box::new(Teno9A),
    };

    let mut writer = FrameWriter::new(cfg.output.clone(), nx, ny)?;
    if writer.cfg.enable {
        // 初期状態の画像も出力しておく
        writer.write_frame(&q, 0.0)?;
    }

    let mut rhs = vec![0.0; nx * ny];
    let mut q1 = vec![0.0; nx * ny];
    let mut q2 = vec![0.0; nx * ny];
    // q3, q4 were used by SSPRK(5,4); not needed for SSPRK(3,3)

    let integrator = cfg.simulation.time_integrator;

    for n in 1..=steps {
        match integrator {
            TimeIntegrator::SspRk3 => {
                // SSPRK(3,3) [1]:
                // q^(1) = q^n + dt * L(q^n)
                scheme_box.rhs(&q, &u, &v, dx, dy, nx, ny, &mut rhs);
                for k in 0..nx * ny {
                    q1[k] = q[k] + dt * rhs[k];
                }

                // q^(2) = 3/4 q^n + 1/4 (q^(1) + dt * L(q^(1)))
                scheme_box.rhs(&q1, &u, &v, dx, dy, nx, ny, &mut rhs);
                for k in 0..nx * ny {
                    q2[k] = 0.75 * q[k] + 0.25 * (q1[k] + dt * rhs[k]);
                }

                // q^{n+1} = 1/3 q^n + 2/3 (q^(2) + dt * L(q^(2)))
                scheme_box.rhs(&q2, &u, &v, dx, dy, nx, ny, &mut rhs);
                for k in 0..nx * ny {
                    q[k] = (1.0 / 3.0) * q[k] + (2.0 / 3.0) * (q2[k] + dt * rhs[k]);
                }
            }
        }

        writer.maybe_write(&q, n, n as f64 * dt)?;
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
    println!(
        "mass0={:.6e} mass1={:.6e} mass_err={:.6e}",
        stats.mass0,
        stats.mass1,
        stats.mass1 - stats.mass0
    );
    println!("L2 = {:.6e}", stats.l2);
    Ok(stats)
}
