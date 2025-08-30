use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use convec::config::{Config, SchemeType, VelocityCfg};
use convec::sim;

fn scan_yaml(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut v = vec![];
    for entry in fs::read_dir(dir).with_context(|| format!("read_dir {}", dir.display()))? {
        let entry = entry?;
        let p = entry.path();
        if p.extension() == Some(OsStr::new("yaml")) {
            v.push(p);
        }
    }
    v.sort();
    Ok(v)
}

fn label(s: SchemeType) -> String {
    format!("{:?}", s).to_lowercase()
}

fn compute_steps(cfg: &Config) -> usize {
    let nx = cfg.simulation.nx;
    let ny = cfg.simulation.ny;
    let lx = cfg.simulation.lx;
    let ly = cfg.simulation.ly;
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;

    // reconstruct velocity and find umax (same as sim::run)
    let mut umax = 0.0f64;
    match cfg.simulation.velocity {
        VelocityCfg::SolidRotation { omega, center_x, center_y } => {
            for j in 0..ny {
                let y = (j as f64 + 0.5) * dy;
                for i in 0..nx {
                    let x = (i as f64 + 0.5) * dx;
                    let u = -omega * (y - center_y);
                    let v = omega * (x - center_x);
                    let s = (u * u + v * v).sqrt();
                    if s > umax { umax = s; }
                }
            }
        }
    }
    let dt = cfg.simulation.cfl * dx.min(dy) / (umax + 1e-15);
    let t_end = cfg.simulation.rotations as f64;
    let steps = (t_end / dt).ceil() as usize;
    // println!("[gallery] steps={} dt={:.3e}", steps, t_end / steps as f64);
    steps
}

fn main() -> Result<()> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let tests_dir = root.join("tests");
    let gallery_dir = root.join("gallery_1rot");
    fs::create_dir_all(&gallery_dir)?;

    let files = scan_yaml(&tests_dir)?;
    for f in files {
        let mut cfg = Config::from_path(&f).with_context(|| format!("load {}", f.display()))?;
        let steps = compute_steps(&cfg);
        let name = label(cfg.scheme.r#type);

        // Prepare a temporary out dir for this run to capture only initial+final
        let tmp_dir = root.join("gallery_tmp").join(name.replace(' ', "_"));
        if tmp_dir.exists() { fs::remove_dir_all(&tmp_dir)?; }
        fs::create_dir_all(&tmp_dir)?;

        // Override output settings
        cfg.output.enable = true;
        cfg.output.dir = tmp_dir.to_string_lossy().to_string();
        cfg.output.prefix = "frame".to_string();
        cfg.output.format = convec::config::OutFmt::Png;
        cfg.output.stride = steps; // write only at final step
        cfg.output.start_index = 0;
        cfg.output.grid = false;
        cfg.output.axes = false;
        cfg.output.colorbar = false;
        // Force colormap to jet for gallery output
        cfg.output.colormap = convec::config::Colormap::Jet;

        // Run
        let _ = sim::run(cfg)?;

        // Find the latest frame in tmp_dir and move as <scheme>.png under gallery_dir
        let mut last_path: Option<PathBuf> = None;
        for entry in fs::read_dir(&tmp_dir)? {
            let p = entry?.path();
            if p.extension() == Some(OsStr::new("png")) {
                last_path = Some(match last_path { Some(cur) => if p > cur { p } else { cur }, None => p });
            }
        }
        if let Some(src) = last_path {
            let dest = gallery_dir.join(format!("{}.png", name));
            fs::rename(&src, &dest).with_context(|| format!("move {} -> {}", src.display(), dest.display()))?;
        }
        // cleanup temp dir
        fs::remove_dir_all(&tmp_dir).ok();
    }
    println!("Saved gallery to {}", gallery_dir.display());
    Ok(())
}
