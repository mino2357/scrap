use std::path::PathBuf;

use convec::{config::Config, sim};

fn run_and_get_l2(file: &str) -> f64 {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(file);
    let cfg = Config::from_path(path).expect("load config");
    let stats = sim::run(cfg).expect("run sim");
    stats.l2
}

#[test]
fn centered8_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/centered8.yaml");
    assert!(l2 < 0.15, "L2 norm too large: {}", l2);
}

#[test]
fn weno5_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/weno5.yaml");
    assert!(l2 < 0.15, "L2 norm too large: {}", l2);
}

#[test]
fn upwind1_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/upwind1.yaml");
    assert!(l2 < 0.6, "L2 norm too large: {}", l2);
}

#[test]
fn tvd_minmod_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/tvd_minmod.yaml");
    assert!(l2 < 0.3, "L2 norm too large: {}", l2);
}

#[test]
fn tvd_vanleer_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/tvd_vanleer.yaml");
    assert!(l2 < 0.25, "L2 norm too large: {}", l2);
}
