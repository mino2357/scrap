set terminal pngcairo size 800,600
set output 'case04.png'
set xlabel 'Time [s]'
set ylabel 'Mole fraction'
set key outside
set grid
set key autotitle columnhead
plot for [i=2:11] 'case04.dat' using 1:i with lines lw 2
