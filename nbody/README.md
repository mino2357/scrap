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
