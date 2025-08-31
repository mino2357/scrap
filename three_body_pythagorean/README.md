# case05: Pythagorean three-body problem

This example integrates Burrau's Pythagorean three-body problem with masses
$3$, $4$, and $5$. Bodies start at rest at the vertices of a $3$-$4$-$5$
right triangle and evolve under their mutual gravitational attraction.

High precision arithmetic using GCC's quad-precision type (`__float128`) and a
custom Bulirsch–Stoer integrator with adaptive time stepping (effective order
$\ge 10$) capture the dynamics accurately. The program writes positions to
`case05.dat` for post-processing.

## Run

```
make run
```

Start gnuplot to display the trajectory:

```
GNUTERM=qt gnuplot plot.gp  # set GNUTERM=dumb for ASCII output in a terminal
```
