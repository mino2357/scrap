# Zalesak DNC (Rust + YAML)

This program solves the two-dimensional passive scalar advection equation

```math
\frac{\partial q}{\partial t} + u\frac{\partial q}{\partial x} + v\frac{\partial q}{\partial y} = 0
```

with a solid body rotation velocity field. Several finite-difference schemes are
available:

* **Centered6** – sixth-order central derivative

  ```math
  \frac{\partial f}{\partial x}\Big|_i \approx \frac{1}{60\Delta x}(f_{i-3}-9 f_{i-2}+45 f_{i-1}-45 f_{i+1}+9 f_{i+2}-f_{i+3})
  ```

* **Centered8** – eighth-order central derivative

  ```math
  \frac{\partial f}{\partial x}\Big|_i \approx \frac{1}{280\Delta x}(-3 f_{i-4}+32 f_{i-3}-168 f_{i-2}+672 f_{i-1}-672 f_{i+1}+168 f_{i+2}-32 f_{i+3}+3 f_{i+4})
  ```

* **Centered10** – tenth-order central derivative
* **Centered12** – twelfth-order central derivative
* **Centered14** – fourteenth-order central derivative
* **WENO5** – fifth-order weighted essentially non-oscillatory scheme
* **WENO5-Z** – improved fifth-order WENO using Z-type weights
* **WENO7-Z** – seventh-order WENO with Z-type weights
* **WENO9-Z** – ninth-order WENO with Z-type weights
* **TENO6** – sixth-order targeted essentially non-oscillatory scheme
* **TENO7-A** – seventh-order TENO variant A
* **TENO8-A** – alias of TENO9-A in this repo (same reconstruction)
* **TENO9-A** – ninth-order TENO variant A
* **Upwind1** – first-order upwind difference
* **Upwind3x3** – first-order upwind difference with a 3×3 stencil

  ```math
  \frac{\partial f}{\partial x}\Big|_{i,j} \approx
  \begin{cases}
    \dfrac{f_{i,j-1}+4 f_{i,j}+f_{i,j+1}-f_{i-1,j-1}-4 f_{i-1,j}-f_{i-1,j+1}}{6\Delta x}, & u_{i,j} \ge 0, \\
    \dfrac{f_{i+1,j-1}+4 f_{i+1,j}+f_{i+1,j+1}-f_{i,j-1}-4 f_{i,j}-f_{i,j+1}}{6\Delta x}, & u_{i,j} < 0.
  \end{cases}
  ```
* **TVD-Minmod**, **TVD-VanLeer** – second-order TVD schemes with flux limiters
* **MP5** – fifth-order monotonicity-preserving scheme

Note: CIP-family schemes (CIP, CIP-CSL, CIP-B, CIP-CSL2, CIP-CSL2-MH) exist in
the codebase but are currently experimental and excluded from the default test
suite and documentation examples.

Time integration uses the third-order TVD Runge–Kutta method of Shu and Osher:

```math
\begin{aligned}
q^{(1)} &= q^n + \Delta t\,L(q^n),\\
q^{(2)} &= \tfrac{3}{4} q^n + \tfrac{1}{4}\big(q^{(1)} + \Delta t\,L(q^{(1)})\big),\\
q^{n+1} &= \tfrac{1}{3} q^n + \tfrac{2}{3}\big(q^{(2)} + \Delta t\,L(q^{(2)})\big).
\end{aligned}
```

Frames are saved with a colour bar, thicker axes, and the current simulation
time overlaid in the corner for easier inspection. The grid overlay is
disabled by default, and the colour map defaults to `jet` (also supporting
`turbo` and `gray`).

## Build & Run
```bash
cargo run --release -- --config config.yaml
```

Switch scheme by editing `scheme.type` (centered6 / centered8 / centered10 /
centered12 / centered14 / weno5 / weno5_z / weno7_z / weno9_z / teno6 /
teno7a / teno8a / teno9a / upwind1 / upwind3x3 / tvd_minmod / tvd_van_leer /
mp5). Frames go to `output.dir`.

## Make video
```bash
ffmpeg -framerate 30 -i frames/cen_%06d.png -pix_fmt yuv420p zalesak.mp4
```
