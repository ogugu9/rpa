set terminal postscript eps enhanced color
set output "eplot.eps"

set size 0.721,1.0
#set size 0.4,1.0

set nokey

#set xrange [-0.2:3.5]
#set yrange [-.2:0.2]
set xtics -1,100,100
#set ytics -4,0.2,4
set xzeroaxis

#plot 'out.eplot.dat','out.eplot.dat' with l 1
plot 'out.eplot.dat' using 4:5 w lp ls 2\
, 'out.eplot.dat' using 1:2 w l ls 1\
#, 'out.eplot.dat' using 1:3 w lp ls 3\

set terminal x11
replot
pause -1
