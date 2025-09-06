# scrap

Numerical experiments for compressible flow and reactive transport.

## nbody: High‑performance N‑Body (CPU/AVX2)

Double‑precision N‑body integrator with AVX2/FMA kernels, Yoshida4/RK4/Verlet, OpenMP, j‑tiling, autosweep, and live plotting via gnuplot.

- Full details live in `nbody/README.md`.

Quick start (from repo root):

- Build: `cd nbody && ./run.sh build --clean`
- Normal run: `./run.sh run --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --method yoshida4 --Bj 512 --threads 4 --kernel fused`
- Live plot: `./run.sh live --N 2048 --steps 2000 --dt 1e-3 --eps 1e-3 --Bj 512 --plot-every 20 --plot-axes xy`
- Compare kernels:
  - Vector vs scalar: `./run.sh compare --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --Bj 512 --threads 4 --kernel fused`
  - Scalar vs scalar_tiled: `./run.sh compare_scalar --N 4096 --steps 50 --dt 1e-3 --eps 1e-3 --Bj 1024 --threads 4 --kernel fused`
- Compiler presets (weak vs strong): `./run.sh compare_opt --N 4096 --steps 200 --dt 1e-3 --eps 1e-3 --Bj 512 --threads 4 --kernel fused`
- Auto‑sweep Bj (scalar_tiled): `./run.sh sweep_scalar_tiled --N 4096 --steps 200 --kernel fused --threads 4 --coarse-min 128 --coarse-max 3072 --coarse-step 128 --fine-radius 256 --fine-step 16 --out sweep_scalar_tiled_4096.csv`

Reading the summary (run.sh run):
- pairs/s = N^2 × steps / elapsed
- GFLOP/s = pairs/s × fpp / 1e9 (default fpp=32)
- GB/s = pairs/s × bytes_per_pair / 1e9 (vector≈8, scalar≈32)
- Arithmetic intensity AI = fpp / bytes_per_pair. See Roofline notes in `nbody/README.md`.

Notes (Windows/WSL2): Foreground apps (e.g. Chrome) can steal budget on low‑power SoCs, reducing pairs/s. Use High‑performance power plan, avoid Efficiency mode on `vmmem`, keep one core free (e.g., threads=Ncores‑1), and prefer static foreground pages during benchmarks. Details in `nbody/README.md`.

## euler1d_sod: 1‑D Euler equations

$\rho$ is density, $u$ velocity, $e$ total energy density and $p$ pressure.
The conservative form of the 1‑D Euler equations is

```math
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x} = 0,
```

```math
\frac{\partial}{\partial t} (\rho u) + \frac{\partial}{\partial x} (p + \rho u^2) = 0,
```

```math
\frac{\partial e}{\partial t} + \frac{\partial}{\partial x} ((e+p)u) = 0,
```

with the equation of state (heat capacity ratio $\gamma$)

```math
p = (\gamma - 1)\left(e - \tfrac{1}{2}\rho u^2\right).
```

The code in `euler1d_sod` reproduces the Sod shock tube problem and compares the
numerical solution against the analytic reference.

## porous_reactor_1d_v1: 1‑D porous catalyst layer

`porous_reactor_1d_v1` solves transport of three species $A,B,C$ and energy in a porous
reactor. Advection with velocity $u$ couples to effective diffusion and a
heterogeneous reaction $A+B\rightarrow C$ on the catalyst surface. The fluid
and solid phases exchange heat and react according to an Arrhenius rate with a
temperature factor $f_T(T_s)$.

Fluid‐phase species balances:

```math
\varepsilon \frac{\partial c_i}{\partial t} + u\frac{\partial c_i}{\partial x}
= D_i\frac{\partial^2 c_i}{\partial x^2} \pm R_{\mathrm{vol}},\quad i=A,B,C,
```

where $R_{\mathrm{vol}}=a_s r_s$ and the surface reaction rate is

```math
r_s = k_0\exp\!\left(-\frac{E_a}{R_g T_s}\right)c_A c_B f_T(T_s).
```

Energy equations for fluid ($T_f$) and solid ($T_s$):

```math
\varepsilon\rho_f C_{p,f}\left(\frac{\partial T_f}{\partial t}+u\frac{\partial T_f}{\partial x}\right)
= \varepsilon k_f\frac{\partial^2 T_f}{\partial x^2} + h_{sf}(T_s-T_f) + \gamma q_{rx},
```

```math
(1-\varepsilon)\rho_s C_{p,s}\frac{\partial T_s}{\partial t}
= (1-\varepsilon)k_s^{\mathrm{eff}}\frac{\partial^2 T_s}{\partial x^2} + h_{sf}(T_f-T_s) + (1-\gamma) q_{rx},
```

with $q_{rx}=a_s(-\Delta H)r_s$. Boundary conditions are Dirichlet at the
inlet and zero‑gradient at the outlet.

## porous_reactor_1d_v2: dimensionless porous reactor

`porous_reactor_1d_v2` nondimensionalises the previous model using length $L$, velocity $u$
and concentration scale $C_0$. Dimensionless groups include

```math
\mathrm{Pe}_i = \frac{uL}{D_i},\quad
\mathrm{Da} = \frac{a_s k_0 C_0 L}{\varepsilon u},\quad
\mathrm{Pe}_{T_f} = \frac{uL}{\alpha_f},\quad
\mathrm{Pe}_{T_s} = \frac{uL}{\alpha_s^{\mathrm{eff}}},
```

along with heat–exchange numbers $H_f$ and $H_s$ and nondimensional reaction
heats $\chi_f$ and $\chi_s$. The dimensionless species and energy equations
mirror those of porous_reactor_1d_v1 but in scaled variables $\hat{x}$, $\hat{t}$, $\hat{c}_i$ and
$\hat{T}_{f,s}$.

## advection2d: 2‑D scalar advection

`advection2d` tracks the rotation of Zalesak's disk under a prescribed velocity
field. Spatial derivatives use either an eighth‑order centred scheme or the
fifth‑order WENO method, while time marching employs a three‑stage TVD
Runge–Kutta integrator. The program writes PNG frames with the current time,
colour bar and labelled axes for later animation. Further details and build
instructions live in `advection2d/README.md`.

## three_body_pythagorean: Pythagorean three-body problem

`three_body_pythagorean` integrates Burrau's Pythagorean three-body problem using GCC's
quad-precision floating point type (`__float128`) and a custom Bulirsch–Stoer
solver with adaptive time stepping (effective order $\ge 10$). Bodies of masses
3, 4, and 5 start at rest at the corners of a 3‑4‑5 triangle and the program
saves their positions in `case05.dat` for later analysis.

## Utilities

The scripts `chemkin_ode/chemkin_to_md.py` and `chemkin_ode/chemkin_to_graph.py` help document CHEMKIN
reaction mechanisms. The former renders the `REACTIONS` block in a Markdown
list, while the latter builds a Graphviz graph of species connectivity. PNG
output from `chemkin_to_graph.py` requires the Graphviz `dot` executable.

## Gaussian rotation benchmark

In addition to Zalesak's slotted disk, a smooth Gaussian bump test is provided
to highlight behaviour on smooth data:

- YAMLs for all schemes live in `advection2d/tests_gaussian/`.
- Print an L2 ranking for one rotation:

  - Zalesak: `cargo run --release --example rank`
  - Gaussian: `cargo run --release --example rank_gauss`

- To render PNG frames for the Gaussian case:

  ```bash
  cargo run --release -- --config advection2d/gaussian.yaml
  ```

### What to expect and why

- On smooth data (Gaussian), high‑order centred schemes are non‑dissipative and
  tend to dominate the L2 ranking: centred14 < centred12 < centred10 < centred8.
- WENO/TENO are designed for robustness near discontinuities; their nonlinear
  weights and upwind fluxes introduce a small but non‑zero dissipation even in
  smooth regions, so their L2 is typically larger for this benchmark.
- Grid/time resolution matters. With a coarse 32×32 grid and 3rd‑order time,
  the advantage of 7th/9th‑order WENO/TENO is reduced. Increasing resolution or
  using the provided 4th‑order SSPRK(5,4) time integrator reveals more of the
  spatial order benefit.

### New options and commands

- WENO/TENO fluxes now use upwind face fluxes: the face velocity
  `u_{i+1/2} = 0.5(u_i+u_{i+1})` picks the upwind reconstruction, which reduces
  dissipation on smooth flows compared to LLF splitting.
- 4th‑order time integrator (SSPRK(5,4)) is available. Enable in YAML:

  ```yaml
  simulation:
    time_integrator: ssp_rk54  # default is ssp_rk3
  ```

- One‑shot galleries for the final frame after one rotation:

  ```bash
  # Zalesak gallery → advection2d/gallery_1rot
  cargo run --release --example gallery

  # Gaussian gallery → advection2d/gallery_gauss_1rot
  cargo run --release --example gallery gauss
  ```

### Notes on TENO variants

- In this repository the current TENO8A/TENO9A implementations reuse the WENO9
  family reconstruction with TENO cut‑off. TENO8A is an alias of TENO9A here,
  so their L2 matches. A distinct 8th‑order coefficient set can be added later
  if desired.

## Comparison of spatial discretisation schemes

The repository includes a variety of one-dimensional advection schemes.  The table
summarises each method's formal accuracy in smooth regions, a ballpark estimate of
per-step cost, and typical uses.  Cost is measured relative to a first-order
upwind method.

| Scheme | Order | Relative cost | Typical use |
|-------|------|---------------|-------------|
| upwind1 | 1 | 1× | cheapest but very diffusive; quick sweeps |
| tvd_minmod | 2 | 1.5× | extremely robust TVD limiter; crude shock capturing |
| tvd_van_leer | 2 | 1.5× | smoother TVD limiter; moderate diffusion |
| upwind3x3 | 3 | 2× | compact 3×3 upwind stencil in 2-D; cheap and stable |
| centered6 | 6 | 2× | low dispersion for smooth waves; no shock control |
| centered8 | 8 | 2.5× | as above with wider stencil |
| centered10 | 10 | 3× | very accurate in periodic problems |
| centered12 | 12 | 3.5× | high-order wave propagation |
| centered14 | 14 | 4× | niche, resolves fine spectra; fragile near shocks |
| mp5 | 5 | 5× | monotonicity preserving; good general-purpose shock capturing |
| weno5 | 5 | 5× | classic WENO-JS; balances cost and accuracy |
| weno5_z | 5 | 5× | WENO-Z weights reduce dissipation |
| weno7_z | 7 | 7× | higher order shock capturing |
| weno9_z | 9 | 9× | very accurate but expensive; research-grade DNS |
| teno6 | 6 | 6× | targeted ENO with selective stencils |
| teno7a | 7 | 7× | adaptive TENO, variant A |
| teno8a | 8 | 8× | alias of TENO9A in this repo |
| teno9a | 9 | 9× | highest-order TENO, variant A |

### Rough ranking

- **Highest fidelity near discontinuities:** TENO9A ≈ TENO8A (alias) ≥ TENO7A ≥ TENO6 ≥ WENO9_Z > WENO7_Z > MP5 ≈ WENO5_Z > WENO5 > TVD_van_leer > TVD_minmod > upwind3x3 > upwind1 > centered schemes.
- **Best for smooth periodic problems:** centered14 ≥ centered12 ≥ centered10 ≥ centered8 ≥ centered6.
- **Cheapest per grid point:** upwind1 < tvd_minmod ≈ tvd_van_leer < upwind3x3 < centered6.

These rankings are intentionally loose but give a flavour of trade-offs among
accuracy, robustness and cost.
