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
