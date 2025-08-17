set terminal pngcairo size 800,600

# Mole fraction history for all species
set output 'case04.png'
set xlabel 'Time [s]'
set ylabel 'Mole fraction'
set logscale y
set key outside
set grid
set key autotitle columnhead
plot for [i=2:11] 'case04.dat' using 1:i with lines lw 2

# Temperature history
set output 'case04_temp.png'
unset logscale y
set ylabel 'Temperature [K]'
plot 'case04.dat' using 1:12 with lines lw 2 title 'T'

# Temperature time derivative
set output 'case04_dTdt.png'
set ylabel 'dT/dt [K/s]'
plot 'case04.dat' using 1:13 with lines lw 2 title 'dT/dt'
