# Zalesak DNC (Rust + YAML)

## Build & Run
```bash
cargo run --release -- --config config.yaml
```

Switch scheme by editing `scheme.type` (centered8 / weno5). Frames go to `output.dir`.

## Make video
```bash
ffmpeg -framerate 30 -i frames/cen_%06d.png -pix_fmt yuv420p zalesak.mp4
```
