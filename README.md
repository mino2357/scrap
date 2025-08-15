# scrap

## case01, 1-dim Euler equations.

$\rho\ \mathrm{(kg/m^3)}$ is fluid (gas) density, $u\ \mathrm{(m/s)}$ is velocity, $e\ \mathrm{(kg/(m \cdot s))}$ is total energy density, $p\ \mathrm{(kg/(m \cdot s^2))}\$ is pressure.

$$
\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x} = 0,
$$

$$
\frac{\partial}{\partial t} \left( \rho u \right) + \frac{\partial}{\partial x} \left( p + \rho u^2 \right) = 0,
$$

$$
\frac{\partial e}{\partial t} + \frac{\partial}{\partial x} \left( \left( e + p \right) u \right) = 0.
$$

Equation of state. $\gamma$ is heat capacity ratio.

$$
p = \left( \gamma - 1 \right) \left( e - \frac{1}{2}\rho u^2 \right).
$$
