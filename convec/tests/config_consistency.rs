use std::ffi::OsStr;
use std::fs;
use std::path::PathBuf;

use convec::config::{Config, InitialConditionCfg, TimeIntegrator, VelocityCfg};

fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[test]
fn all_yaml_share_common_zalesak_settings() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let dir = root.join("tests");

    let mut checked = 0usize;
    for entry in fs::read_dir(&dir).expect("read tests dir") {
        let entry = entry.expect("dir entry");
        let path = entry.path();
        if path.extension() != Some(OsStr::new("yaml")) {
            continue;
        }

        let name = path.file_name().unwrap().to_string_lossy().to_string();
        let cfg = Config::from_path(&path).expect("load yaml config");

        // Simulation grid and CFL
        assert_eq!(cfg.simulation.nx, 32, "{}: nx", name);
        assert_eq!(cfg.simulation.ny, 32, "{}: ny", name);
        assert!(approx(cfg.simulation.lx, 1.0, 1e-12), "{}: lx", name);
        assert!(approx(cfg.simulation.ly, 1.0, 1e-12), "{}: ly", name);
        assert!(approx(cfg.simulation.cfl, 0.40, 1e-12), "{}: cfl", name);
        assert_eq!(cfg.simulation.rotations, 1, "{}: rotations", name);

        // Time integrator: default SSPRK(3,3) with one exception allowed
        // - centered8_ssprk54.yaml uses SSPRK(5,4)
        let expected_ti = if name == "centered8_ssprk54.yaml" {
            TimeIntegrator::SspRk54
        } else {
            TimeIntegrator::SspRk3
        };
        assert!(
            matches!((cfg.simulation.time_integrator, expected_ti),
                (TimeIntegrator::SspRk3, TimeIntegrator::SspRk3)
                | (TimeIntegrator::SspRk54, TimeIntegrator::SspRk54)
            ),
            "{}: unexpected time integrator: {:?}",
            name,
            cfg.simulation.time_integrator
        );

        // Velocity: solid rotation with omega=2π centered at (0.5, 0.5)
        match cfg.simulation.velocity {
            VelocityCfg::SolidRotation { omega, center_x, center_y } => {
                let two_pi = std::f64::consts::PI * 2.0;
                assert!(approx(omega, two_pi, 1e-12), "{}: omega", name);
                assert!(approx(center_x, 0.5, 1e-12), "{}: vel.center_x", name);
                assert!(approx(center_y, 0.5, 1e-12), "{}: vel.center_y", name);
            }
        }

        // Initial condition: Zalesak disk with fixed parameters
        match cfg.initial_condition {
            InitialConditionCfg::Zalesak { center_x, center_y, radius, slot_width, slot_length } => {
                assert!(approx(center_x, 0.5, 1e-12), "{}: ic.center_x", name);
                assert!(approx(center_y, 0.75, 1e-12), "{}: ic.center_y", name);
                assert!(approx(radius, 0.18, 1e-12), "{}: ic.radius", name);
                assert!(approx(slot_width, 0.06, 1e-12), "{}: ic.slot_width", name);
                assert!(approx(slot_length, 0.26, 1e-12), "{}: ic.slot_length", name);
            }
            _ => panic!("{}: expected Zalesak initial condition", name),
        }

        checked += 1;
    }

    // Sanity: ensure we actually checked a decent number of YAMLs
    assert!(checked >= 10, "expected to check many YAML test cases, got {}", checked);
}
