use std::ffi::OsStr;
use std::fs;
use std::path::PathBuf;

use convec::{config::{Config, VelocityCfg}, sim};

fn run_and_get_l2_with_sign(path: &PathBuf, sign: f64) -> f64 {
    let mut cfg = Config::from_path(path).expect("load config");
    // Flip rotation direction by changing omega's sign.
    cfg.simulation.velocity = match cfg.simulation.velocity {
        VelocityCfg::SolidRotation { omega, center_x, center_y } => VelocityCfg::SolidRotation {
            omega: sign * omega.abs(),
            center_x,
            center_y,
        },
    };
    let stats = sim::run(cfg).expect("run sim");
    stats.l2
}

fn scan_yaml(dir: &PathBuf) -> Vec<PathBuf> {
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

// Use a modest relative tolerance to account for minor asymmetries
// in upwinding and time-stepping on coarse grids.
const TOL_REL: f64 = 0.05; // 5%

#[test]
#[ignore]
fn gaussian_reverse_rotation_l2_similar_across_schemes() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let dir = root.join("tests_gaussian");
    let files = scan_yaml(&dir);
    assert!(!files.is_empty(), "no gaussian YAMLs found in {:?}", dir);

    for yaml in files {
        let l2_fwd = run_and_get_l2_with_sign(&yaml, 1.0);
        let l2_rev = run_and_get_l2_with_sign(&yaml, -1.0);
        let denom = l2_fwd.max(l2_rev).max(1e-12);
        let rel = (l2_fwd - l2_rev).abs() / denom;
        assert!(
            rel <= TOL_REL,
            "L2 mismatch for {}: forward={:.6e}, reverse={:.6e}, rel_diff={:.3}% (tol {:.1}%)",
            yaml.display(),
            l2_fwd,
            l2_rev,
            rel * 100.0,
            TOL_REL * 100.0
        );
    }
}

#[test]
fn gaussian_centered6_reverse_rotation_l2_similar() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let yaml = root.join("tests_gaussian/centered6.yaml");
    let l2_fwd = run_and_get_l2_with_sign(&yaml, 1.0);
    let l2_rev = run_and_get_l2_with_sign(&yaml, -1.0);
    let denom = l2_fwd.max(l2_rev).max(1e-12);
    let rel = (l2_fwd - l2_rev).abs() / denom;
    assert!(
        rel <= TOL_REL,
        "L2 mismatch for {}: forward={:.6e}, reverse={:.6e}, rel_diff={:.3}% (tol {:.1}%)",
        yaml.display(),
        l2_fwd,
        l2_rev,
        rel * 100.0,
        TOL_REL * 100.0
    );
}

#[test]
fn gaussian_weno5z_reverse_rotation_l2_similar() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let yaml = root.join("tests_gaussian/weno5_z.yaml");
    let l2_fwd = run_and_get_l2_with_sign(&yaml, 1.0);
    let l2_rev = run_and_get_l2_with_sign(&yaml, -1.0);
    let denom = l2_fwd.max(l2_rev).max(1e-12);
    let rel = (l2_fwd - l2_rev).abs() / denom;
    assert!(
        rel <= TOL_REL,
        "L2 mismatch for {}: forward={:.6e}, reverse={:.6e}, rel_diff={:.3}% (tol {:.1}%)",
        yaml.display(),
        l2_fwd,
        l2_rev,
        rel * 100.0,
        TOL_REL * 100.0
    );
}

#[test]
fn gaussian_teno8a_reverse_rotation_l2_similar() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let yaml = root.join("tests_gaussian/teno8a.yaml");
    let l2_fwd = run_and_get_l2_with_sign(&yaml, 1.0);
    let l2_rev = run_and_get_l2_with_sign(&yaml, -1.0);
    let denom = l2_fwd.max(l2_rev).max(1e-12);
    let rel = (l2_fwd - l2_rev).abs() / denom;
    assert!(
        rel <= TOL_REL,
        "L2 mismatch for {}: forward={:.6e}, reverse={:.6e}, rel_diff={:.3}% (tol {:.1}%)",
        yaml.display(),
        l2_fwd,
        l2_rev,
        rel * 100.0,
        TOL_REL * 100.0
    );
}
