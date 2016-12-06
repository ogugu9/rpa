set terminal aqua dashed
set colorsequence podo
set size 0.7,1
set nokey
set xrange [0: 1.75]
set yrange [ -3 :  3]
set arrow from  0.000,  -3 to  0.000,   3 nohead
set arrow from  0.250,  -3 to  0.250,   3 nohead
set arrow from  0.500,  -3 to  0.500,   3 nohead
set arrow from  1.000,  -3 to  1.000,   3 nohead
set arrow from  1.250,  -3 to  1.250,   3 nohead
set arrow from  1.750,  -3 to  1.750,   3 nohead
set arrow from  0.000,   0.000 to  1.750,   0.000 nohead
set xtics (" G "  0.000," X "  0.250," M "  0.500," G "  1.00," X "  1.250," X "  1.750) font "Times New Roman,22"
set ytics ("{/=20 {/Times-Italic E}^F}" 0) font "Times New Roman,20"
set ytics 1
set ylabel "{/=20 {/Times-Italic E} (eV)}"
p 'out.band.dn.dat' every :1 u 1:2 w l lw 1
