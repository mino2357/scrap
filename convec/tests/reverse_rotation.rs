use std::path::PathBuf;

use convec::{config::{Config, VelocityCfg}, sim};

fn run_and_get_l2_with_sign(file: &str, sign: f64) -> f64 {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(file);
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

fn assert_l2_almost_equal_forward_reverse(yaml: &str, tol_rel: f64) {
    let l2_fwd = run_and_get_l2_with_sign(yaml, 1.0);
    let l2_rev = run_and_get_l2_with_sign(yaml, -1.0);
    let denom = l2_fwd.max(l2_rev).max(1e-12);
    let rel = (l2_fwd - l2_rev).abs() / denom;
    assert!(
        rel <= tol_rel,
        "L2 mismatch for {}: forward={:.6e}, reverse={:.6e}, rel_diff={:.3}% (tol {:.1}%)",
        yaml,
        l2_fwd,
        l2_rev,
        rel * 100.0,
        tol_rel * 100.0
    );
}

// Use a modest relative tolerance to account for minor asymmetries
// in upwinding and time-stepping on coarse grids.
const TOL_REL: f64 = 0.05; // 5%

#[test]
fn centered6_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/centered6.yaml", TOL_REL);
}

#[test]
fn centered8_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/centered8.yaml", TOL_REL);
}

#[test]
fn centered10_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/centered10.yaml", TOL_REL);
}

#[test]
fn centered12_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/centered12.yaml", TOL_REL);
}

#[test]
fn centered14_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/centered14.yaml", TOL_REL);
}

#[test]
fn weno5_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/weno5.yaml", TOL_REL);
}

#[test]
fn weno5_z_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/weno5_z.yaml", TOL_REL);
}

#[test]
fn weno7_z_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/weno7_z.yaml", TOL_REL);
}

#[test]
fn weno9_z_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/weno9_z.yaml", TOL_REL);
}

#[test]
fn teno6_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/teno6.yaml", TOL_REL);
}

#[test]
fn teno7a_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/teno7a.yaml", TOL_REL);
}

#[test]
fn teno8a_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/teno8a.yaml", TOL_REL);
}

#[test]
fn teno9a_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/teno9a.yaml", TOL_REL);
}

#[test]
fn upwind1_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/upwind1.yaml", TOL_REL);
}

#[test]
fn upwind3x3_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/upwind3x3.yaml", TOL_REL);
}

#[test]
fn tvd_minmod_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/tvd_minmod.yaml", TOL_REL);
}

#[test]
fn tvd_vanleer_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/tvd_vanleer.yaml", TOL_REL);
}

#[test]
fn mp5_reverse_rotation_l2_similar() {
    assert_l2_almost_equal_forward_reverse("tests/mp5.yaml", TOL_REL);
}

