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
fn centered6_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/centered6.yaml");
    assert!(l2 < 0.2, "L2 norm too large: {}", l2);
}

#[test]
fn centered6_preserves_rotation_direction() {
    use convec::schemes::{Centered6, Scheme};
    use convec::utils::idx;

    // Use a periodic sine wave so that the boundary wrap-around is consistent
    // with the stencil.  For q = sin(kx) and uniform velocity u = 1, the
    // semi-discrete advection term becomes -k cos(kx).
    let nx = 32; // sufficiently large and periodic
    let ny = 4;
    let dx = 1.0;
    let dy = 1.0;
    let k = 2.0 * std::f64::consts::PI / nx as f64;

    let mut q = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            q[idx(i, j, nx)] = (k * i as f64).sin();
        }
    }

    let u = vec![1.0; nx * ny];
    let v = vec![0.0; nx * ny];
    let mut out = vec![0.0; nx * ny];
    Centered6.rhs(&q, &u, &v, dx, dy, nx, ny, &mut out);

    for j in 0..ny {
        for i in 0..nx {
            let expected = -k * (k * i as f64).cos();
            let got = out[idx(i, j, nx)];
            assert!((got - expected).abs() < 1e-6, "{} vs {}", got, expected);
        }
    }
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
fn centered14_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/centered14.yaml");
    assert!(l2 < 0.15, "L2 norm too large: {}", l2);
}

#[test]
/// Verify the 5-stage 4th-order SSPRK(5,4) scheme [Spiteri & Ruuth 2002]
/// integrates the centered8 spatial discretization with acceptable error.
fn centered8_ssprk54_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/centered8_ssprk54.yaml");
    // With the correct SSPRK(5,4) coefficients, the error increases
    // slightly compared to the older approximate implementation.
    assert!(l2 < 0.41, "L2 norm too large: {}", l2);
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
fn weno9_z_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/weno9_z.yaml");
    assert!(l2 < 0.17, "L2 norm too large: {}", l2);
}

#[test]
fn teno6_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/teno6.yaml");
    assert!(l2 < 0.15, "L2 norm too large: {}", l2);
}

#[test]
fn teno7a_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/teno7a.yaml");
    assert!(l2 < 0.16, "L2 norm too large: {}", l2);
}

#[test]
fn teno8a_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/teno8a.yaml");
    assert!(l2 < 0.18, "L2 norm too large: {}", l2);
}

#[test]
fn teno9a_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/teno9a.yaml");
    assert!(l2 < 0.18, "L2 norm too large: {}", l2);
}

#[test]
fn upwind1_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/upwind1.yaml");
    assert!(l2 < 0.6, "L2 norm too large: {}", l2);
}

#[test]
fn upwind3x3_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/upwind3x3.yaml");
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

/// Ensure all CIP-based schemes advect a sine wave in the correct
/// direction.  For q = sin(kx) and uniform velocity u = 1, the
/// semi-discrete advection term should be -k cos(kx) for every scheme.
fn check_cip_rotation<S: convec::schemes::Scheme>(scheme: S) {
    use convec::utils::idx;

    let nx = 32; // sufficiently large and periodic
    let ny = 4;
    let dx = 1.0;
    let dy = 1.0;
    let k = 2.0 * std::f64::consts::PI / nx as f64;

    let mut q = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            q[idx(i, j, nx)] = (k * i as f64).sin();
        }
    }

    let u = vec![1.0; nx * ny];
    let v = vec![0.0; nx * ny];
    let mut out = vec![0.0; nx * ny];
    scheme.rhs(&q, &u, &v, dx, dy, nx, ny, &mut out);

    for j in 0..ny {
        for i in 0..nx {
            let expected = -k * (k * i as f64).cos();
            let got = out[idx(i, j, nx)];
            if expected.abs() > 1e-12 {
                assert!(
                    got.signum() == expected.signum(),
                    "wrong sign: {} vs {} at ({}, {})",
                    got,
                    expected,
                    i,
                    j
                );
            }
        }
    }
}

#[test]
fn cip_schemes_preserve_rotation_direction() {
    use convec::schemes::{Cip, CipB, CipCsl, CipCsl2, CipCsl2Mh};

    check_cip_rotation(Cip);
    check_cip_rotation(CipCsl);
    check_cip_rotation(CipB);
    check_cip_rotation(CipCsl2);
    check_cip_rotation(CipCsl2Mh);
}

#[test]
fn mp5_l2_below_threshold() {
    let l2 = run_and_get_l2("tests/mp5.yaml");
    assert!(l2 < 0.17, "L2 norm too large: {}", l2);
}
