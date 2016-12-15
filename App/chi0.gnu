set terminal aqua dashed
set colorsequence podo
set size 0.7,1
set nokey
set xrange [0: 3.701]
set yrange [ -0.02 :  0.04]
set arrow from  0.000,  -0.02 to  0.000,   0.04 nohead
set arrow from  0.102,  -0.02 to  0.102,   0.04 nohead
set arrow from  1.049,  -0.02 to  1.049,   0.04 nohead
set arrow from  1.152,  -0.02 to  1.152,   0.04 nohead
set arrow from  1.275,  -0.02 to  1.275,   0.04 nohead
set arrow from  2.083,  -0.02 to  2.083,   0.04 nohead
set arrow from  2.186,  -0.02 to  2.186,   0.04 nohead
set arrow from  2.682,  -0.02 to  2.682,   0.04 nohead
set arrow from  0.000,   0.000 to  3.701,   0.000 nohead
set xtics (" G "  0.000," Z "  0.102," B "  1.049," D "  1.152," S "  1.275," Y "  2.083," T "  2.186," R "  2.682," G "  3.701) font "Times New Roman,22"
set ytics ("{/=20 {/Times-Italic E}^F}" 0) font "Times New Roman,20"
set ytics 0.01
set ylabel "{/=20 {/Times-Italic E} (eV)}"
p 'chi0.dat' u 1:2 w l lw 1
