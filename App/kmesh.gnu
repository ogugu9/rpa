set terminal aqua dashed
set colorsequence podo
set size 1,1
set nokey
set xrange [ -1 :  1]
set yrange [ -1 :  1]
set arrow from  0.000,  -4.000 to  0.000,   4.000 nohead
set arrow from -4.000,   0.000 to  4.000,   0.000 nohead
p 'kmesh.dat' u 5:6, 'sympt.dat' u 1:2 w p pt 1 ps 2 lc 7
