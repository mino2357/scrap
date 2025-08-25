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
fn centered10_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/centered10.yaml");
    assert!(l2 < 0.15, "L2 norm too large: {}", l2);
}

#[test]
fn centered12_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/centered12.yaml");
    assert!(l2 < 0.15, "L2 norm too large: {}", l2);
}

#[test]
/// Verify the 5-stage 4th-order SSPRK(5,4) scheme [Spiteri & Ruuth 2002]
/// integrates the centered8 spatial discretization with acceptable error.
fn centered8_ssprk54_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/centered8_ssprk54.yaml");
    assert!(l2 < 0.35, "L2 norm too large: {}", l2);
}

#[test]
fn weno5_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/weno5.yaml");
    assert!(l2 < 0.15, "L2 norm too large: {}", l2);
}

#[test]
fn weno5_z_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/weno5_z.yaml");
    assert!(l2 < 0.13, "L2 norm too large: {}", l2);
}

#[test]
fn weno7_z_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/weno7_z.yaml");
    assert!(l2 < 0.17, "L2 norm too large: {}", l2);
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

#[test]
fn cip_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/cip.yaml");
    assert!(l2 < 0.3, "L2 norm too large: {}", l2);
}

#[test]
fn cip_csl_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/cip_csl.yaml");
    assert!(l2 < 0.3, "L2 norm too large: {}", l2);
}

#[test]
fn cip_b_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/cip_b.yaml");
    assert!(l2 < 0.3, "L2 norm too large: {}", l2);
}

#[test]
fn cip_csl2_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/cip_csl2.yaml");
    assert!(l2 < 0.31, "L2 norm too large: {}", l2);
}

#[test]
fn cip_csl2_mh_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/cip_csl2_mh.yaml");
    assert!(l2 < 0.31, "L2 norm too large: {}", l2);
}

#[test]
fn mp5_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/mp5.yaml");
    assert!(l2 < 0.17, "L2 norm too large: {}", l2);
}
