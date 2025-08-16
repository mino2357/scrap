# Sample CHEMKIN-II Solver

This example reads a small CHEMKIN-II mechanism (`chem.inp`) and integrates the
resulting set of ordinary differential equations.  Reaction rates follow
mass–action kinetics and the temporal evolution of the species is obtained via an
embedded Runge–Kutta (RK45) method with adaptive step sizing.

### Initial conditions

The simulation assumes a constant temperature of **1000 K** and integrates from
`t = 0` to `t = 1e-3 s`.  Species are initialized in moles as

| species | value |
|---------|-------|
| H2      | 1.0   |
| H       | 0.5   |
| others  | 1e-8  |

The output file normalizes these values to mole fractions so that the sum of all
species at each time is one.

### Visualization

Running the example produces `case04.dat` containing time histories of mole
fractions.  A bundled gnuplot script (`plot.gp`) is invoked automatically to
create `case04.png`, visualizing mole fraction versus time for every species.

## Build and run

```bash
make -C case04 run
```

The command builds the solver, runs the simulation and generates the PNG plot.
