use anyhow::{Context, Result};
use std::env;
use std::path::{Path, PathBuf};

use convec::config::Config;
use convec::sim;

/// エントリポイント。
///
/// `--config <path>` で設定ファイルを指定可能。相対パスは
/// 実行ディレクトリに加えて `CARGO_MANIFEST_DIR`（クレートルート）
/// 相対も探索する。
/// 設定ファイルの探索順序:
/// 1) 与えられたパスをそのまま
/// 2) クレートルートからの相対
/// 3) `crate_name/` プレフィクスを剥がしてクレートルート相対
fn resolve_cfg_path(p: &Path) -> PathBuf {
    // Try as-is first
    if p.is_absolute() || p.exists() {
        return p.to_path_buf();
    }
    // Try relative to crate root (CARGO_MANIFEST_DIR)
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cand = manifest_dir.join(p);
    if cand.exists() {
        return cand;
    }
    // If path is prefixed by crate dir name (e.g. "convec/gaussian.yaml")
    if let Some(crate_name) = manifest_dir.file_name() {
        if let Ok(stripped) = p.strip_prefix(crate_name) {
            let cand2 = manifest_dir.join(stripped);
            if cand2.exists() {
                return cand2;
            }
        }
    }
    // Fallback: return original
    p.to_path_buf()
}

fn main() -> Result<()> {
    // デフォルトの設定ファイルは同ディレクトリの `config.yaml`
    let mut cfg_path = PathBuf::from("config.yaml");

    // `env::args` でプログラム名を除いた引数を取得し，
    // `--config <path>` もしくは `-c <path>` が指定されていれば
    // そのパスを設定ファイルとして用いる。
    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        if a == "--config" || a == "-c" {
            if let Some(p) = args.next() {
                cfg_path = PathBuf::from(p);
            }
        }
    }

    // 設定ファイルを読み込み，各種情報を表示する。
    // Resolve common invocation patterns (workspace root or crate dir)
    let cfg_path = resolve_cfg_path(&cfg_path);
    println!("[info] loading config: {}", cfg_path.display());
    let cfg = Config::from_path(&cfg_path).with_context(|| format!("failed to read {}", cfg_path.display()))?;
    println!("[info] {:?}", cfg.summary());

    // 実際のシミュレーションを実行。戻り値は現在未使用だが
    // 将来的に統計情報等へ利用できる。
    let _stats = sim::run(cfg)?;
    Ok(())
}
