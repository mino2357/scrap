use anyhow::{Context, Result};
use serde::Deserialize;
use std::fs;
use std::path::Path;

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case")]
pub struct Config {
    pub simulation: SimulationCfg,
    pub initial_condition: InitialConditionCfg,
    pub scheme: SchemeCfg,
    pub output: OutputCfg,
}

impl Config {
    pub fn from_path<P: AsRef<Path>>(p: P) -> Result<Self> {
        let s = fs::read_to_string(&p)
            .with_context(|| format!("failed to read {}", p.as_ref().display()))?;
        let cfg: Config = serde_yaml::from_str(&s).with_context(|| "YAML parse error")?;
        Ok(cfg)
    }
    pub fn summary(&self) -> String {
        format!(
            "scheme={:?}, time_int={:?}, N=({},{}) CFL={} ROT={} out_dir={} stride={} fmt={:?}",
            self.scheme.r#type,
            self.simulation.time_integrator,
            self.simulation.nx,
            self.simulation.ny,
            self.simulation.cfl,
            self.simulation.rotations,
            self.output.dir,
            self.output.stride,
            self.output.format
        )
    }
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case")]
pub struct SimulationCfg {
    pub nx: usize,
    pub ny: usize,
    pub lx: f64,
    pub ly: f64,
    pub cfl: f64,
    pub rotations: usize,
    pub velocity: VelocityCfg,
    #[serde(default)]
    pub time_integrator: TimeIntegrator,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum VelocityCfg {
    SolidRotation {
        omega: f64,
        center_x: f64,
        center_y: f64,
    },
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum InitialConditionCfg {
    Zalesak {
        center_x: f64,
        center_y: f64,
        radius: f64,
        slot_width: f64,
        slot_length: f64,
    },
    Disk {
        center_x: f64,
        center_y: f64,
        radius: f64,
    },
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case")]
pub struct SchemeCfg {
    pub r#type: SchemeType,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum SchemeType {
    Centered8,
    Centered10,
    Centered12,
    Weno5,
    Upwind1,
    Cip,
    CipCsl,
    CipB,
    CipCsl2,
    CipCsl2Mh,
    Mp5,
    TvdMinmod,
    TvdVanLeer,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
/// Strong-stability-preserving Runge–Kutta time integrators.
/// - `SspRk3`: 3-stage 3rd-order scheme of Shu & Osher (1989).
/// - `SspRk54`: 5-stage 4th-order scheme of Spiteri & Ruuth (2002).
pub enum TimeIntegrator {
    /// SSPRK(3,3)
    SspRk3,
    /// SSPRK(5,4)
    SspRk54,
}

impl Default for TimeIntegrator {
    fn default() -> Self {
        TimeIntegrator::SspRk3
    }
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case")]
pub struct OutputCfg {
    pub enable: bool,
    pub dir: String,
    pub prefix: String,
    pub format: OutFmt,
    pub stride: usize,
    pub start_index: usize,
    pub scale: ScaleCfg,
    pub flip_y: bool,
    pub out_w: Option<usize>,
    pub out_h: Option<usize>,
    pub interp: Interp,
    pub grid: bool,
    pub grid_step: usize,
    pub grid_thick: usize,
    pub axes: bool,
    pub colormap: Colormap,
    pub colorbar: bool,
}

#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case", tag = "mode")]
pub enum ScaleCfg {
    Fixed { min: f64, max: f64 },
    Auto,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum OutFmt {
    Png,
    Ppm,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum Interp {
    Nearest,
    Bilinear,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum Colormap {
    Gray,
    Turbo,
}
