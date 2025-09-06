# N-Body (DP) — Yoshida 4th-Order Symplectic vs RK4 (WSL2/Windows)

This is a minimal, **double-precision** N-body code for CPUs with **AVX2+FMA** (e.g., Intel N97).  
It includes:
- Force kernel: AVX2 (DP) with **SP rsqrt + DP Newton×2** (for fast `1/sqrt(r^2)` on AVX2).
- **Yoshida 4th-order symplectic (1990)** integrator (DKD), **RK4**, and **Velocity-Verlet (2nd)**.
- SoA layout, L1-friendly **j-tiling**, OpenMP parallelism.
- Built‑in **self check** (`--selftest`) running a 2-body circular-orbit sanity test.
- Energy & angular momentum diagnostics.

> Works great under **WSL2**; also compiles on native Windows with **MSVC** (or clang‑cl).

---

## Build (WSL2 / Linux, Clang recommended)

```bash
sudo apt update
sudo apt install -y clang lld cmake ninja-build build-essential libomp-dev
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_EXE_LINKER_FLAGS="-fuse-ld=lld"
cmake --build build -j
```

## Build (Windows, MSVC / clang-cl)

- Install **Visual Studio 2022** (Desktop C++), CMake, Ninja.
- Open *x64 Native Tools Command Prompt for VS 2022*.
```bat
cmake -S . -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```
MSVC flags are set in `CMakeLists.txt` (`/O2 /arch:AVX2 /fp:fast`). For clang‑cl the generator will use similar flags.

---

## Run

```bash
# Example: 4096 bodies, Yoshida4, 200 steps, dt=1e-3, softening eps=1e-3, j-tile=512
./nbody --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512 --threads 4 --check 1
```

Key options:
- `--N <int>`: number of bodies (any positive integer).
- `--steps <int>`: time steps.
- `--dt <float>`: time step size (constant; required for symplectic behavior).
- `--eps <float>`: Plummer softening (epsilon).
- `--method {yoshida4|rk4|verlet}`.
- `--Bj <int>`: j‑tile size (default 512; tune to L1).
- `--threads <int>`: OpenMP threads (N97 = 4).
- `--check {0|1}`: compute energy/angmom periodically.
- `--seed <uint64>`: RNG seed.
- `--selftest`: runs a 2‑body circular‑orbit test and prints relative errors.
- `--kernel {fused|two_pass}`: fused kick or two-pass forces.
- `--mode {vector|scalar|scalar_tiled}`:
  - `vector`: AVX2/FMA + SP rsqrt + NR×2 + jタイル（既定）
  - `scalar`: 素朴なスカラ（二重ループ, タイル無し）
  - `scalar_tiled`: スカラ＋j方向タイル（Bj）だけ適用（SIMDなし）

A typical Yoshida4 run prints pairs/s, GFLOP/s (est.), and relative energy & angular‑momentum errors.

---

## Helper script (`run.sh`)

You can use the unified helper to build and run, and to print a detailed throughput summary at the end:

```bash
# Build (once)
./run.sh build --clean

# Normal run (same as calling ./nbody ...), plus a summary with AI/GB/s
./run.sh run --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 \
             --method yoshida4 --Bj 512 --threads 4 --kernel fused --check 1 \
             --bytes-per-pair 8 --fpp 32

# Optional: show roofline context by providing machine peaks
./run.sh run ... --peak-gflops 115 --mem-gbs 40

# Live plotting via gnuplot
./run.sh live --N 2048 --steps 2000 --dt 1e-3 --eps 1e-3 --Bj 512 --plot-every 20
```

Flags `--bytes-per-pair` (default 8) and `--fpp` (default 32 FLOP/pair) control the summary’s memory bandwidth estimate and GFLOP/s conversion. Provide `--peak-gflops` and `--mem-gbs` to print roofline limits and utilization.

### Compare compiler optimization presets

Build presets are selectable via CMake cache `NBODY_OPT_PRESET` and `run.sh build --opt`:

- default: `-O3 -march=native -mavx2 -mfma -ffp-contract=fast -fopenmp`
- strong: `-Ofast -march=native -mavx2 -mfma -ffast-math -ffp-contract=fast -fopenmp -flto -fomit-frame-pointer -funroll-loops`
- weak: `-O0 -g -march=native -mavx2 -mfma -fopenmp -fno-builtin -fno-unroll-loops -fno-omit-frame-pointer`

Quick A/B with separate build dirs and speedup report:

```bash
./run.sh compare_opt --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 \
                     --method yoshida4 --Bj 512 --threads 4 --kernel fused
```

Manual control of a single build dir:

```bash
# Strong
./run.sh build --opt strong
# Weak
./run.sh build --opt weak
```

### Compare naive vs optimized kernels

Run the same case twice (vector vs scalar) and report speedup:

```bash
./run.sh compare --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 \
                 --method yoshida4 --Bj 512 --threads 4 --kernel fused
```

Tip for summary consistency:
- Vector mode: `--bytes-per-pair 8` (j=32B shared by i=4 lanes)
- Scalar mode: `--bytes-per-pair 32`（スカラは共有がないため）

### Compare scalar vs scalar_tiled (タイル効果の分離)

スカラ実装に Bj タイルだけを入れた `scalar_tiled` と素朴スカラ `scalar` を比較して、SIMD ではなくアルゴリズム（局所性改善）の寄与だけを見る:

```bash
./run.sh compare_scalar --N 4096 --steps 50 --dt 1e-3 --eps 1e-3 \
                        --method yoshida4 --Bj 1024 --threads 4 --kernel fused
```

注意: この差はケース依存で小さいことがあります。N と Bj によっては HW プリフェッチが効き、`scalar` が十分に速く見えることもあります。Bj を掃引（例: 512, 1024, 1536, 2048）し、`N=8192` 以上など L2/DRAM を跨ぐ領域で観察すると違いが見えやすくなります。集計時の `--bytes-per-pair` はスカラ系では 32 を目安にしてください。

### scalar_tiled の最適 Bj を自動探索（N97想定）

`run.sh sweep` は Bj の粗→細スイープを行い、最良値を表示します。`MODE=scalar_tiled` を付けるとスカラ＋タイルの探索ができます。

例（4096体, 200ステップ, Bj=128..3072 粗→±256を16刻みで微調整）:

```bash
# 長時間の測定向け（コマンドは nbody ディレクトリで実行）
MODE=scalar_tiled N=4096 STEPS=200 COARSE_MIN=128 COARSE_MAX=3072 \
COARSE_STEP=128 FINE_RADIUS=256 FINE_STEP=16 OUT=sweep_scalar_tiled_4096.csv \
./run.sh sweep --kernel fused --threads 4
```

ヒント（N97/Gracemontのメモリ事情を踏まえて）:
- j 粒子のワーキングセットは 1 粒子 32B（x,y,z,m）。Bj の塊あたり `32×Bj` Bytes。
- 「タイルを変えた直後の i が直前の j を再利用できる」ことが効くので、L1D に“直近の j 塊”が残る範囲が狙い目です。
- 実測では Bj≈1.5k〜2.0k 付近で“台地”が出やすい一方、scalar_tiled は lane 共有がないため plateu が緩く、Bj の感度がケース依存になります。粗→細の自動探索で最良を拾ってください。
- スループットの読み方は本 README の「Throughput summary」を参照。スカラ系は `--bytes-per-pair 32` を基準にメモリ帯域の目安を読むと、Roofline の位置づけがしやすいです。

---

## Reading the throughput summary (formulas and rationale)

When you run via `./run.sh run ...`, a summary like the following is printed:

```text
[SUMMARY]
- config: method=yoshida4 kernel=fused Bj=512 threads=4 eps=0.001
- size: N=4096 steps=200
- timing: elapsed=10.417s  pairs/s=3.221000e+08
- pairs: total=3.355443e+09
- compute: fpp=32  GFLOP/s=10.307
- memory: bytes/pair=8  est GB/s=2.577  AI=4.00 flop/byte
```

Meaning of each field and how to compute it:

- config: method, kernel, `Bj`, `threads`, `eps` used in the run.
- size: problem size parameters `N` and `steps`.
- timing: `elapsed` is wall time in seconds; `pairs/s` is the interaction rate.
  - Formula: `pairs/s = (N^2 × steps) / elapsed`.
- pairs: `total = N^2 × steps` is the total number of i–j interactions executed.
- compute: FLOP rate using an assumed FLOP per pair (fpp).
  - Formula: `GFLOP/s = (pairs/s × fpp) / 1e9`.
  - Default `fpp = 32` is a common rough budget for a softened gravitational pair: r², ε² add, reciprocal square root, inv³, and 3× FMA accumulations; exact counts vary by kernel.
- memory: bandwidth estimate using effective bytes per pair (bpp) and arithmetic intensity (AI).
  - Formulas:
    - `GB/s = (pairs/s × bpp) / 1e9`
    - `AI = fpp / bpp` (FLOP per byte)
  - Default `bpp = 8` models the current AVX2 fused kernel where each j‑particle load is 32 bytes (xj,yj,zj,mj as doubles) broadcast/shared across i=4 lanes, i.e. `32 B / 4 = 8 B` per i–j pair. If you disable lane sharing (scalar demo), set `--bytes-per-pair 32`.

Worked example for the snippet above (`N=4096`, `steps=200`, `elapsed=10.417 s`, `pairs/s=3.221×10^8`, `fpp=32`, `bpp=8`):

```math
\text{pairs}_\text{total} = N^2\,\times\,\text{steps} = 4096^2\times 200 = 3.3554432\times10^9
```

```math
\text{pairs/s} = \frac{\text{pairs}_\text{total}}{\text{elapsed}} \approx \frac{3.3554\times10^9}{10.417} \approx 3.221\times10^8
```

```math
\text{GFLOP/s} = \frac{\text{pairs/s}\times fpp}{10^9} = \frac{3.221\times10^8\times 32}{10^9} \approx 10.31
```

```math
\text{GB/s} = \frac{\text{pairs/s}\times bpp}{10^9} = \frac{3.221\times10^8\times 8}{10^9} \approx 2.58
\quad\text{and}\quad
AI = \frac{fpp}{bpp} = \frac{32}{8} = 4.0\;\text{FLOP/byte}
```

Optional: roofline context. If you pass `--peak-gflops` (compute peak) and `--mem-gbs` (sustainable memory bandwidth), the runner prints roofline limits and utilization. The model is

```math
P_{\max} = \min\bigl(\,\pi,\; \beta \times AI\,\bigr),
```

where `π` is the peak compute throughput (GFLOP/s), `β` the memory bandwidth (GB/s), and `AI` the arithmetic intensity. The achieved GFLOP/s compared to `π` gives compute‑peak utilization; the estimated GB/s compared to `β` gives memory‑bandwidth utilization.

References
- Roofline model: S. Williams, A. Waterman, D. Patterson, “Roofline: An Insightful Visual Performance Model for Multicore Architectures” (Commun. ACM 52(4), 2009; LBNL Tech Report). `P = min(π, β×I)`.
- Instruction throughput/latency (for understanding rsqrt/sqrt cost, FMAs, etc.): Agner Fog, “Instruction tables: Lists of latencies, throughputs and micro‑operation breakdowns”.
- Yoshida 4th‑order symplectic integrator: H. Yoshida, “Construction of higher order symplectic integrators” (Physics Letters A 150, 1990).

---

## Notes
- Arrays are **64B aligned**; AVX2 kernel uses `_mm_rsqrt_ps` for SP initial guess + DP Newton×2.
- Self‑interaction is benign (r=0 ⇒ alpha·r = 0), so not explicitly skipped.
- FTZ/DAZ enabled at start to avoid denorm slowdowns.
- For best WSL2 performance: keep repo under Linux FS (e.g., `~/`), pin threads (`OMP_PLACES=cores OMP_PROC_BIND=close`).

### Windows/WSL2: Foreground app hurts throughput (why and mitigations)

Symptom: When switching focus to a foreground GUI app (e.g. Chrome), pairs/s can drop noticeably.

Why this happens on low‑power SoCs (e.g. Intel N97):
- Foreground boost: Windows gives the foreground app a scheduling boost; background work (WSL’s `vmmem`) may get shorter quanta and lower priority.
- Power throttling/EcoQoS: Background apps can be shifted to energy‑efficient QoS, which lowers CPU frequency aggressiveness (Windows 11 “Efficiency mode”).
- Shared power budget: iGPU activity (browser rendering/video) and CPU share a tight power envelope; GPU load can pull down CPU turbo clocks.
- WSL2 VM boundary: Linux nice/affinity applies inside the VM; Windows ultimately schedules `vmmem`, which is treated as a background process if the terminal/VSCode isn’t focused.

Mitigations (pick what fits your setup):
- Power plan: Set Windows power mode to “Best performance” and select the “High performance” or “Ultimate Performance” plan. In Advanced settings → Processor power management, set Minimum processor state=100%, Cooling policy=Active. Reboot or `wsl --shutdown` after changing.
- Task Manager: Ensure `vmmem`/`VmmemWSL` is not in “Efficiency mode”. Optionally set its Priority to High (Details → right‑click → Set priority → High).
- Reserve a core for UI: run with one fewer thread (e.g. `OMP_NUM_THREADS=3` on a 4‑core N97) so the OS always has headroom for the foreground app.
- Keep the browser idle: avoid animated/video tabs during benchmarks; keep a blank tab foreground. Hardware acceleration ON usually reduces CPU load; if GPU saturation drags CPU clocks, prefer a static page.
- Reduce GUI in WSL: for long runs avoid live plotting (`--plot none`), or increase `--plot_every`.
- After toggling settings: `wsl --shutdown` to restart the VM so new host policies apply cleanly.

References: Windows “Efficiency mode” (Task Manager) and EcoQoS; Microsoft’s Foreground Boost behavior; Roofline reasoning above explains the impact of CPU/GPU power sharing on compute‑bound kernels.

---

## License
MIT


## Live plotting (gnuplot)
Install gnuplot inside WSL2:
```bash
sudo apt install -y gnuplot
```
Run with live view (updates every 20 steps, plotting XY plane, up to 4096 points):
```bash
./nbody --N 2048 --steps 2000 --dt 1e-3 --eps 1e-3 \
        --method yoshida4 --Bj 512 --threads 4 \
        --plot live --plot_every 20 --plot_limit 4096 --plot_axes xy
```
Make sure WSLg (GUI) is active so the gnuplot window appears.
