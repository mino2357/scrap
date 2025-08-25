use anyhow::Result;
use std::env;
use std::path::PathBuf;

use convec::config::Config;
use convec::sim;

/// エントリポイント。
///
/// コマンドライン引数から設定ファイルのパスを読み取り，
/// シミュレーションを実行する。
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
    println!("[info] loading config: {}", cfg_path.display());
    let cfg = Config::from_path(&cfg_path)?;
    println!("[info] {:?}", cfg.summary());

    // 実際のシミュレーションを実行。戻り値は現在未使用だが
    // 将来的に統計情報等へ利用できる。
    let _stats = sim::run(cfg)?;
    Ok(())
}
