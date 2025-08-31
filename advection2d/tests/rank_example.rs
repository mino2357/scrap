use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};

use convec::config::{Config, SchemeType};
use convec::sim;

fn scan_yaml(dir: &Path) -> Vec<PathBuf> {
    let mut v = vec![];
    for entry in fs::read_dir(dir).expect("read_dir") {
        let entry = entry.expect("dir entry");
        let p = entry.path();
        if p.extension() == Some(OsStr::new("yaml")) {
            v.push(p);
        }
    }
    v.sort();
    v
}

fn label(s: SchemeType) -> String {
    format!("{:?}", s).to_lowercase()
}

// Full rotation ranking across all YAMLs.
#[test]
#[ignore]
fn rotation_ranking_full() {
    // Use the same directory as the example to gather cases
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let dir = root.join("tests");
    let files = scan_yaml(&dir);

    // Run all and collect L2 by scheme label
    let mut results: Vec<(String, f64)> = vec![];
    for f in files {
        let cfg = Config::from_path(&f).expect("load config");
        let stats = sim::run(cfg.clone()).expect("run sim");
        let name = label(cfg.scheme.r#type);
        assert!(stats.l2.is_finite() && stats.l2 >= 0.0, "non-finite L2 for {}", name);
        results.push((name, stats.l2));
    }
    // Sort ascending and ensure non-decreasing order
    results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    for w in results.windows(2) {
        assert!(w[0].1 <= w[1].1 + 1e-15, "ranking not non-decreasing");
    }

    // Check alias behaviour: teno8a and teno9a produce identical L2
    let mut l2_8a = None;
    let mut l2_9a = None;
    for (name, l2) in &results {
        if name == "teno8a" {
            l2_8a = Some(*l2);
        } else if name == "teno9a" {
            l2_9a = Some(*l2);
        }
    }
    let l2_8a = l2_8a.expect("missing teno8a case");
    let l2_9a = l2_9a.expect("missing teno9a case");
    assert!((l2_8a - l2_9a).abs() <= 1e-15, "teno8a vs teno9a mismatch: {} vs {}", l2_8a, l2_9a);

    // A couple of loose ordering sanity checks that reflect the example ranking
    let mut map = std::collections::HashMap::new();
    for (name, l2) in results {
        map.insert(name, l2);
    }
    if let (Some(w9z), Some(tvdl)) = (map.get("weno9z"), map.get("tvdvanleer")) {
        assert!(w9z < tvdl, "expected weno9_z to beat tvd_van_leer: {} vs {}", w9z, tvdl);
    }
    if let (Some(tvdl), Some(tvdm)) = (map.get("tvdvanleer"), map.get("tvdminmod")) {
        assert!(tvdl < tvdm, "expected tvd_van_leer to beat tvd_minmod: {} vs {}", tvdl, tvdm);
    }
}
