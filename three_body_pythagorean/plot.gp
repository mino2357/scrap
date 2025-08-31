set title 'Pythagorean three-body problem'
set xlabel 'x'
set ylabel 'y'
set size ratio -1
set key left top
do for [i=0:1e9] {
    plot 'case05.dat' using 2:3 with lines lc rgb 'red' title 'm1', \
         '' using 4:5 with lines lc rgb 'green' title 'm2', \
         '' using 6:7 with lines lc rgb 'blue' title 'm3'
    pause 0.1
}
