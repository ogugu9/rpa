set terminal postscript eps enhanced color
set output 'fv.eps'

set nokey 
set size square 

#set xrange[-0.5:0.5]
#set yrange[-0.5:0.5]
set xrange[-0.7:0.7]
set yrange[-0.7:0.7]
set xtics -0.5,0.25,0.5
set ytics -0.5,0.25,0.5

#set xtics -0.5,0.025,0.5
#set ytics -0.5,0.025,0.5

set style arrow 1 nohead lw 1
set arrow as 1

set arrow 1 from -0.5,-0.5 to (0.5),(-0.5) nohead back lw 0
set arrow 2 from (0.5),(-0.5) to (0.5),(0.5) nohead back lw 0
set arrow 3 from (0.5),(0.5) to (-0.5),(0.5) nohead back lw 0
set arrow 4 from (-0.5),(0.5) to (-0.5),(-0.5) nohead back lw 0

a= 0.1

plot 'v-out.fs.dat'using 1:2:($3)*a:($4)*a with vector,\
'out.fs.dat'using 5:6:($7-$5):($8-$6) with vector
# as 1


#set grid 
#plot 'out.fs.dat'using 5:6:7:8 every ::1::5 w vector as 1
#plot 'gomi.dat'w vector as 1

set output
set term x11
#set term wxt
replot 
pause -1

exit

