set terminal postscript eps enhanced color
set output 'fermi.eps'

set key outside
set size square

set xrange[-0.5:0.5]
set yrange[-0.5:0.5]
set xtics -0.5,0.25,0.5
set ytics -0.5,0.25,0.5

plot 'out.fermi.dat' using 1:2 w d notitle "1",\
'out.fermi.dat' using 3:4 w p title "2",\
'out.fermi.dat' using 5:6 w p title "3",\
'out.fermi.dat' using 7:8 w p title "4",\
'out.fermi.dat' using 9:10 w p title "5",\
'out.fermi.dat' using 11:12 w p title "6",\
'out.fermi.dat' using 13:14 w p title "7",\
'out.fermi.dat' using 15:16 w p title "8",\
'out.fermi.dat' using 17:18 w p title "9",\
'out.fermi.dat' using 19:20 w p title "10",\
'out.fermi.dat' using 21:22 w p title "11",\
'out.fermi.dat' using 23:24 w p title "12",\
'out.fermi.dat' using 25:26 w p title "13",\
'out.fermi.dat' using 27:28 w p title "14 ",\
'out.fermi.dat' using 29:30 w p title "15",\
'out.fermi.dat' using 31:32 w p title "16",\
'out.fermi.dat' using 33:34 w p title "17",\
'out.fermi.dat' using 35:36 w p title "18",\
'out.fermi.dat' using 37:38 w p title "19",\
'out.fermi.dat' using 39:40 w p title "20"


set output
set term x11
#set term wxt
replot 
pause -1

exit




set pm3d map

splot 'out.fermi2.dat'

set output
set term x11
replot 
pause -1
