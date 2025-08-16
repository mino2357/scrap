# Sample CHEMKIN-II Solver

This example reads a small CHEMKIN-II mechanism (`chem.inp`) and integrates the resulting
set of ordinary differential equations using an embedded Runge-Kutta (RK45) method with
adaptive step sizing.

## Build and run

```bash
make -C case04 run
```
