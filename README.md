# scrap

Numerical experiments for compressible flow and reactive transport.

## case01: 1‑D Euler equations

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

The code in `case01` reproduces the Sod shock tube problem and compares the
numerical solution against the analytic reference.

## case02: 1‑D porous catalyst layer

`case02` solves transport of three species $A,B,C$ and energy in a porous
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

## case03: dimensionless porous reactor

`case03` nondimensionalises the previous model using length $L$, velocity $u$
and concentration scale $C_0$. Dimensionless groups include

```math
\mathrm{Pe}_i = \frac{uL}{D_i},\quad
\mathrm{Da} = \frac{a_s k_0 C_0 L}{\varepsilon u},\quad
\mathrm{Pe}_{T_f} = \frac{uL}{\alpha_f},\quad
\mathrm{Pe}_{T_s} = \frac{uL}{\alpha_s^{\mathrm{eff}}},
```

along with heat–exchange numbers $H_f$ and $H_s$ and nondimensional reaction
heats $\chi_f$ and $\chi_s$. The dimensionless species and energy equations
mirror those of case02 but in scaled variables $\hat{x}$, $\hat{t}$, $\hat{c}_i$ and
$\hat{T}_{f,s}$.

## convec: 2‑D scalar advection

`convec` tracks the rotation of Zalesak's disk under a prescribed velocity
field. Spatial derivatives use either an eighth‑order centred scheme or the
fifth‑order WENO method, while time marching employs a three‑stage TVD
Runge–Kutta integrator. The program writes PNG frames with the current time,
colour bar and labelled axes for later animation. Further details and build
instructions live in `convec/README.md`.

## case05: Pythagorean three-body problem

`case05` integrates Burrau's Pythagorean three-body problem using GCC's
quad-precision floating point type (`__float128`) and a custom Bulirsch–Stoer
solver with adaptive time stepping (effective order $\ge 10$). Bodies of masses
3, 4, and 5 start at rest at the corners of a 3‑4‑5 triangle and the program
saves their positions in `case05.dat` for later analysis.

## Utilities

The scripts `chemkin_to_md.py` and `chemkin_to_graph.py` help document CHEMKIN
reaction mechanisms. The former renders the `REACTIONS` block in a Markdown
list, while the latter builds a Graphviz graph of species connectivity. PNG
output from `chemkin_to_graph.py` requires the Graphviz `dot` executable.
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
| cip | 3 | 3× | cubic profile; needs stored gradients |
| cip_csl | 4 | 3.5× | conservative CIP using cubic spline |
| cip_b | 4 | 4× | B-spline flavour; smoother at extra cost |
| cip_csl2 | 4 | 4× | extended conservative scheme |
| cip_csl2_mh | 4 | 5× | multi-moment hybrid; most accurate CIP variant |
| mp5 | 5 | 5× | monotonicity preserving; good general-purpose shock capturing |
| weno5 | 5 | 5× | classic WENO-JS; balances cost and accuracy |
| weno5_z | 5 | 5× | WENO-Z weights reduce dissipation |
| weno7_z | 7 | 7× | higher order shock capturing |
| weno9_z | 9 | 9× | very accurate but expensive; research-grade DNS |
| teno6 | 6 | 6× | targeted ENO with selective stencils |
| teno7a | 7 | 7× | adaptive TENO, variant A |
| teno8a | 8 | 8× | higher-order TENO, variant A |
| teno9a | 9 | 9× | highest-order TENO, variant A |

### Rough ranking

- **Highest fidelity near discontinuities:** TENO9A ≥ TENO8A ≥ TENO7A ≥ TENO6 ≥ WENO9_Z > WENO7_Z > MP5 ≈ WENO5_Z > WENO5 > CIP_CSL2_MH > CIP_CSL2 > CIP_B ≈ CIP_CSL > CIP > TVD_van_leer > TVD_minmod > upwind3x3 > upwind1 > centered schemes.
- **Best for smooth periodic problems:** centered14 ≥ centered12 ≥ centered10 ≥ centered8 ≥ centered6.
- **Cheapest per grid point:** upwind1 < tvd_minmod ≈ tvd_van_leer < upwind3x3 < centered6.

These rankings are intentionally loose but give a flavour of trade-offs among
accuracy, robustness and cost.
