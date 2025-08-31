use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use convec::config::{Config, SchemeType};
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

fn main() -> Result<()> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let dir = root.join("tests");
    let files = scan_yaml(&dir)?;
    let mut results: Vec<(String, f64)> = vec![];
    for f in files {
        let cfg = Config::from_path(&f).with_context(|| format!("load {}", f.display()))?;
        // Run and collect L2
        let stats = sim::run(cfg.clone()).with_context(|| "run sim")?;
        let name = label(cfg.scheme.r#type);
        results.push((name, stats.l2));
    }
    results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    println!("# L2 ranking (ascending) for Zalesak, 1 rotation");
    for (i, (name, l2)) in results.iter().enumerate() {
        println!("{:2}. {:<18} L2={:.6e}", i + 1, name, l2);
    }
    Ok(())
}
