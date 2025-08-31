set datafile separator ','
set terminal pngcairo size 800,600 enhanced
set output 'sod.png'
set title 'Sod shock tube at t=0.2'
set xlabel 'x'
set ylabel 'Density'
plot 'solution.csv' using 1:2 with lines title 'numerical', \
     'exact.csv' using 1:2 with lines title 'exact'
