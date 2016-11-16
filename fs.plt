set terminal postscript eps enhanced color
set output 'fs.eps'

set key outside
set size square

set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set xtics -0.5,0.25,0.5
set ytics -0.5,0.25,0.5

set style arrow 1 nohead lw 1
set arrow as 1

plot 'out.fs.dat'using 5:6:($7-$5):($8-$6) with vector as 1

set output
set term x11
#set term wxt
replot 
pause -1

exit

