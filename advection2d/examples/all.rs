use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use convec::config::{Colormap, Config, OutFmt, SchemeType, VelocityCfg};
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

fn flip_omega(cfg: &mut Config) {
    cfg.simulation.velocity = match cfg.simulation.velocity {
        VelocityCfg::SolidRotation { omega, center_x, center_y } => VelocityCfg::SolidRotation {
            omega: -omega.abs(),
            center_x,
            center_y,
        },
    };
}

fn rank(dir: &Path, reverse: bool, title: &str) -> Result<()> {
    let files = scan_yaml(dir)?;
    let mut results: Vec<(String, f64)> = vec![];
    for f in files {
        let mut cfg = Config::from_path(&f).with_context(|| format!("load {}", f.display()))?;
        if reverse { flip_omega(&mut cfg); }
        let stats = sim::run(cfg.clone()).with_context(|| "run sim")?;
        let name = label(cfg.scheme.r#type);
        results.push((name, stats.l2));
    }
    results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    println!("{}", title);
    for (i, (name, l2)) in results.iter().enumerate() {
        println!("{:2}. {:<18} L2={:.6e}", i + 1, name, l2);
    }
    println!("");
    Ok(())
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
    steps
}

fn gallery(dir: &Path, out_dir: &Path, reverse: bool) -> Result<()> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    fs::create_dir_all(out_dir)?;
    let files = scan_yaml(dir)?;
    for f in files {
        let mut cfg = Config::from_path(&f).with_context(|| format!("load {}", f.display()))?;
        let steps = compute_steps(&cfg);
        let name = label(cfg.scheme.r#type);

        let tmp_dir = root.join("gallery_tmp").join(name.replace(' ', "_"));
        if tmp_dir.exists() { fs::remove_dir_all(&tmp_dir)?; }
        fs::create_dir_all(&tmp_dir)?;

        cfg.output.enable = true;
        cfg.output.dir = tmp_dir.to_string_lossy().to_string();
        cfg.output.prefix = "frame".to_string();
        cfg.output.format = OutFmt::Png;
        cfg.output.stride = steps;
        cfg.output.start_index = 0;
        cfg.output.grid = false;
        cfg.output.axes = false;
        cfg.output.colorbar = false;
        cfg.output.colormap = Colormap::Jet;
        if reverse { flip_omega(&mut cfg); }

        let _ = sim::run(cfg)?;

        let mut last_path: Option<PathBuf> = None;
        for entry in fs::read_dir(&tmp_dir)? {
            let p = entry?.path();
            if p.extension() == Some(OsStr::new("png")) {
                last_path = Some(match last_path { Some(cur) => if p > cur { p } else { cur }, None => p });
            }
        }
        if let Some(src) = last_path {
            let dest = out_dir.join(format!("{}.png", name));
            fs::rename(&src, &dest).with_context(|| format!("move {} -> {}", src.display(), dest.display()))?;
        }
        fs::remove_dir_all(&tmp_dir).ok();
    }
    Ok(())
}

fn main() -> Result<()> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let tests = root.join("tests");
    let tests_gauss = root.join("tests_gaussian");
    let gal_z_fwd = root.join("gallery_1rot");
    let gal_z_rev = root.join("gallery_1rot_rev");
    let gal_g_fwd = root.join("gallery_gauss_1rot");
    let gal_g_rev = root.join("gallery_gauss_1rot_rev");

    // Rankings
    rank(&tests, false, "# L2 ranking (Zalesak, 1 rot)")?;
    rank(&tests, true,  "# L2 ranking (Zalesak, 1 rot, REVERSED)")?;
    rank(&tests_gauss, false, "# L2 ranking (Gaussian, 1 rot)")?;
    rank(&tests_gauss, true,  "# L2 ranking (Gaussian, 1 rot, REVERSED)")?;

    // Galleries
    gallery(&tests, &gal_z_fwd, false)?;
    gallery(&tests, &gal_z_rev, true)?;
    gallery(&tests_gauss, &gal_g_fwd, false)?;
    gallery(&tests_gauss, &gal_g_rev, true)?;

    println!("Saved galleries to:\n  {}\n  {}\n  {}\n  {}", gal_z_fwd.display(), gal_z_rev.display(), gal_g_fwd.display(), gal_g_rev.display());
    Ok(())
}

