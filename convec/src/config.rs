use anyhow::{Context, Result};
use serde::Deserialize;
use std::fs;
use std::path::Path;

/// シミュレーション全体の設定を保持するトップレベル構造体。
///
/// `config.yaml` から読み込まれ，各モジュールへ設定値を渡す。
#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case")]
pub struct Config {
    /// 計算格子や CFL 数など時間発展に関する設定
    pub simulation: SimulationCfg,
    /// 初期条件の形状などを定義
    pub initial_condition: InitialConditionCfg,
    /// 空間離散化スキームの選択
    pub scheme: SchemeCfg,
    /// 画像出力に関する設定
    pub output: OutputCfg,
}

impl Config {
    /// YAML ファイルから `Config` を生成するユーティリティ。
    pub fn from_path<P: AsRef<Path>>(p: P) -> Result<Self> {
        let s = fs::read_to_string(&p)
            .with_context(|| format!("failed to read {}", p.as_ref().display()))?;
        let cfg: Config = serde_yaml::from_str(&s).with_context(|| "YAML parse error")?;
        Ok(cfg)
    }
    /// ログ出力用の簡易サマリを返す。
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

/// 格子数や計算領域，時間積分方法を含む数値シミュレーションに関する設定。
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

/// 速度場の設定。現在は剛体回転のみをサポートする。
#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum VelocityCfg {
    /// 原点または任意の点を中心とした剛体回転
    SolidRotation {
        omega: f64,
        center_x: f64,
        center_y: f64,
    },
}

/// 初期スカラー場の形状を表す設定。
#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case", tag = "type")]
pub enum InitialConditionCfg {
    /// Zalesak's disk
    Zalesak {
        center_x: f64,
        center_y: f64,
        radius: f64,
        slot_width: f64,
        slot_length: f64,
    },
    /// 単純な円盤
    Disk {
        center_x: f64,
        center_y: f64,
        radius: f64,
    },
}

/// 空間差分スキームの選択。
#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case")]
pub struct SchemeCfg {
    pub r#type: SchemeType,
}

/// 利用可能な空間差分スキーム。
#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum SchemeType {
    Centered6,
    Centered8,
    Centered10,
    Centered12,
    Centered14,
    Weno5,
    Weno5Z,
    Weno7Z,
    Weno9Z,
    Upwind1,
    Upwind3x3,
    Cip,
    CipCsl,
    CipB,
    CipCsl2,
    CipCsl2Mh,
    Mp5,
    TvdMinmod,
    TvdVanLeer,
    Teno6,
    #[serde(rename = "teno7a")]
    Teno7A,
    #[serde(rename = "teno8a")]
    Teno8A,
    #[serde(rename = "teno9a")]
    Teno9A,
}

#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
/// Strong-stability-preserving Runge–Kutta time integrators.
/// - `SspRk3`: 3-stage 3rd-order scheme of Shu & Osher (1989).
pub enum TimeIntegrator {
    /// SSPRK(3,3)
    SspRk3,
}

impl Default for TimeIntegrator {
    fn default() -> Self {
        TimeIntegrator::SspRk3
    }
}

/// 画像としての結果出力に関する設定。
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
    #[serde(default)]
    pub grid: bool,
    #[serde(default = "default_grid_step")]
    pub grid_step: usize,
    #[serde(default = "default_grid_thick")]
    pub grid_thick: usize,
    pub axes: bool,
    #[serde(default)]
    pub colormap: Colormap,
    #[serde(default)]
    pub colorbar: bool,
}

/// カラーマップのスケーリング方法。
#[derive(Debug, Deserialize, Clone)]
#[serde(rename_all = "snake_case", tag = "mode")]
pub enum ScaleCfg {
    /// 最小値と最大値を固定
    Fixed { min: f64, max: f64 },
    /// データに合わせて自動スケーリング
    Auto,
}

/// 出力画像のフォーマット。
#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum OutFmt {
    Png,
    Ppm,
}

/// 出力画像生成時の補間方法。
#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum Interp {
    Nearest,
    Bilinear,
}

/// 利用可能なカラーマップ。
#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(rename_all = "snake_case")]
pub enum Colormap {
    Gray,
    Turbo,
    Jet,
}

impl Default for Colormap {
    fn default() -> Self {
        Colormap::Jet
    }
}

fn default_grid_step() -> usize {
    8
}

fn default_grid_thick() -> usize {
    1
}
