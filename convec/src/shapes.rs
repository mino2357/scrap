use crate::config::InitialConditionCfg;
use crate::utils::idx;

/// 初期条件を生成する。
///
/// Zalesak のスロット付き円柱テスト\[1\]は，移流スキームの形状保存性を評価する
/// 代表的なベンチマークとして広く利用されている。
/// 本関数ではその形状を含むいくつかの初期条件を生成する。
///
/// [1] S. T. Zalesak, "Fully multidimensional flux-corrected transport algorithms for fluids",
/// *Journal of Computational Physics*, 31(3), 335-362, 1979.
pub fn init_field(cfg: &InitialConditionCfg, nx: usize, ny: usize, lx: f64, ly: f64) -> Vec<f64> {
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let mut q = vec![0.0; nx * ny];
    match cfg {
        InitialConditionCfg::Zalesak {
            center_x,
            center_y,
            radius,
            slot_width,
            slot_length,
        } => {
            for j in 0..ny {
                let y = (j as f64 + 0.5) * dy;
                for i in 0..nx {
                    let x = (i as f64 + 0.5) * dx;
                    let r2 = (x - center_x) * (x - center_x) + (y - center_y) * (y - center_y);
                    let disk = r2 <= radius * radius;
                    let slot = ((x - center_x).abs() <= slot_width * 0.5)
                        && (y < center_y + radius)
                        && (y > center_y + radius - slot_length);
                    q[idx(i, j, nx)] = if disk && !slot { 1.0 } else { 0.0 };
                }
            }
        }
        InitialConditionCfg::Disk {
            center_x,
            center_y,
            radius,
        } => {
            for j in 0..ny {
                let y = (j as f64 + 0.5) * dy;
                for i in 0..nx {
                    let x = (i as f64 + 0.5) * dx;
                    let r2 = (x - center_x) * (x - center_x) + (y - center_y) * (y - center_y);
                    q[idx(i, j, nx)] = if r2 <= radius * radius { 1.0 } else { 0.0 };
                }
            }
        }
        InitialConditionCfg::Gaussian {
            center_x,
            center_y,
            sigma,
            amplitude,
        } => {
            let sig2 = sigma * sigma;
            for j in 0..ny {
                let y = (j as f64 + 0.5) * dy;
                for i in 0..nx {
                    let x = (i as f64 + 0.5) * dx;
                    let r2 = (x - center_x) * (x - center_x) + (y - center_y) * (y - center_y);
                    q[idx(i, j, nx)] = amplitude * (-r2 / (2.0 * sig2)).exp();
                }
            }
        }
    }
    q
}
