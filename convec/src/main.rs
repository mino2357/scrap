use anyhow::Result;
use std::env;
use std::path::PathBuf;

mod config;
mod utils;
mod shapes;
mod render;
mod schemes;
mod sim;

use config::Config;

fn main() -> Result<()> {
    let mut cfg_path = PathBuf::from("config.yaml");
    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        if a == "--config" || a == "-c" {
            if let Some(p) = args.next() { cfg_path = PathBuf::from(p); }
        }
    }
    println!("[info] loading config: {}", cfg_path.display());
    let cfg = Config::from_path(&cfg_path)?;
    println!("[info] {:?}", cfg.summary());
    sim::run(cfg)?;
    Ok(())
}
