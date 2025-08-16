# Sample CHEMKIN-II Solver

This example reads a small CHEMKIN-II mechanism (`chem.inp`) and integrates the
resulting set of ordinary differential equations.  Reaction rates follow
mass–action kinetics and the temporal evolution of the species and temperature
is obtained via explicit Runge–Kutta schemes.  The calculation assumes a
constant pressure of **202650 Pa**.  The program supports a classical
fourth‑order method as well as adaptive Runge–Kutta–Fehlberg 4(5) and 7(8)
integrators.

### Initial conditions

The simulation integrates from `t = 0` to `t = 1e-3 s` starting at **1000 K**
and **202650 Pa**.  Species are initialized (in arbitrary moles) as

| species | value |
|---------|-------|
| H2      | 1.0   |
| O2      | 1.0   |
| N2      | 3.76  |
| others  | 1e-8  |

The solver normalizes these values to mole fractions, computes concentrations
from the ideal‑gas law at the specified pressure and temperature, and evolves
species and temperature using NASA polynomial thermodynamic data loaded from
`therm.dat`.  The output file (`case04.dat`) lists mole fractions of all
species and the temperature history.

### Visualization

Running the example produces `case04.dat` containing time histories of mole
fractions and temperature.  A bundled gnuplot script (`plot.gp`) is invoked
automatically to create `case04.png`, visualizing mole fraction versus time for
every species.

## Build and run

```bash
make -C case04 run
```

The command builds the solver, runs the simulation and generates the PNG plot.
To experiment with other integrators, invoke the executable directly, e.g.

```bash
./case04/chem rk45   # Runge–Kutta–Fehlberg 4(5)
./case04/chem rk78   # Runge–Kutta–Fehlberg 7(8)
```
