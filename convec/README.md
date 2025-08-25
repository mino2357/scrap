# Zalesak DNC (Rust + YAML)

This program solves the two-dimensional passive scalar advection equation

```math
\frac{\partial q}{\partial t} + u\frac{\partial q}{\partial x} + v\frac{\partial q}{\partial y} = 0
```

with a solid body rotation velocity field. Several finite-difference schemes are
available:

* **Centered8** – eighth-order central derivative
* **Centered10** – tenth-order central derivative
* **Centered12** – twelfth-order central derivative
* **Centered2d** – second-order central derivative using a 3×3 stencil

  ```math
  \frac{\partial f}{\partial x}\Big|_i \approx \frac{1}{280\Delta x}(-3 f_{i-4}+32 f_{i-3}-168 f_{i-2}+672 f_{i-1}-672 f_{i+1}+168 f_{i+2}-32 f_{i+3}+3 f_{i+4})
  ```

* **WENO5** – fifth-order weighted essentially non-oscillatory scheme
* **WENO5-Z** – improved fifth-order WENO using Z-type weights
* **WENO7-Z** – seventh-order WENO with Z-type weights
* **Upwind1** – first-order upwind difference
* **Upwind2d** – first-order upwind difference with a 3×3 stencil
* **TVD-Minmod**, **TVD-VanLeer** – second-order TVD schemes with flux limiters
* **CIP**, **CIP-CSL**, **CIP-B** – cubic interpolated propagation schemes

  ```math
  \omega_k = \frac{\alpha_k}{\sum_{m=0}^2 \alpha_m},\qquad \alpha_k = \frac{d_k}{(\varepsilon+\beta_k)^2}
  ```

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

Switch scheme by editing `scheme.type` (centered8 / centered10 / centered12 / weno5 / upwind1 / tvd_minmod / tvd_van_leer). Frames go to `output.dir`.

## Make video
```bash
ffmpeg -framerate 30 -i frames/cen_%06d.png -pix_fmt yuv420p zalesak.mp4
```
