set terminal aqua dashed
set colorsequence podo
set size 0.7,1
set nokey
set xrange [0: 3.701]
set yrange [ -2 :  4]
set arrow from  0.000,  -2 to  0.000,   4 nohead
set arrow from  0.102,  -2 to  0.102,   4 nohead
set arrow from  1.049,  -2 to  1.049,   4 nohead
set arrow from  1.152,  -2 to  1.152,   4 nohead
set arrow from  1.275,  -2 to  1.275,   4 nohead
set arrow from  2.083,  -2 to  2.083,   4 nohead
set arrow from  2.186,  -2 to  2.186,   4 nohead
set arrow from  2.682,  -2 to  2.682,   4 nohead
set arrow from  0.000,   0.000 to  3.701,   0.000 nohead
set xtics (" G "  0.000," Z "  0.102," B "  1.049," D "  1.152," S "  1.275," Y "  2.083," T "  2.186," R "  2.682," G "  3.701) font "Times New Roman,22"
set ytics ("{/=20 {/Times-Italic E}^F}" 0) font "Times New Roman,20"
set ylabel "{/=20 {/Times-Italic E} (eV)}"
p 'out.band.dn.dat' every :1 u 1:2 w l lw 1
