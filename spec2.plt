
#set size square
set size 1,0.5

set xrange[-2.:2.5]
set xrange[-0.5:0.5]
#set xrange[-0.5:0.5]
#set yrange[-0.5:0.5]
set zeroaxis



set terminal postscript eps enhanced color
set output 'spec2.eps'
set title "down spin on (1,1)"

plot 'out.spec2.dat' using 1:2 w l title "3z^2-r^2",\
     'out.spec2.dat' using 1:3 w l title "zx",\
     'out.spec2.dat' using 1:4 w l title "yz",\
     'out.spec2.dat' using 1:5 w l title "x^2-y^2",\
     'out.spec2.dat' using 1:6 w l title "xy"


set term x11
set output
#set term wxt
replot 
pause -1

set terminal postscript eps enhanced color
set output 'spec2.eps'
set title "down spin on (2,1)"

plot 'out.spec2.dat' using 1:7  w l title "3z^2-r^2",\
     'out.spec2.dat' using 1:8  w l title "zx",\
     'out.spec2.dat' using 1:9  w l title "yz",\
     'out.spec2.dat' using 1:10 w l title "x^2-y^2",\
     'out.spec2.dat' using 1:11 w l title "xy"            



set term x11
set output
#set term wxt
replot 
pause -1


set terminal postscript eps enhanced color
set output 'spec2.eps'
set title "down spin on (3,1)"

plot 'out.spec2.dat' using 1:12  w l title "3z^2-r^2",\
     'out.spec2.dat' using 1:13  w l title "zx",\
     'out.spec2.dat' using 1:14  w l title "yz",\
     'out.spec2.dat' using 1:15 w l title "x^2-y^2",\
     'out.spec2.dat' using 1:16 w l title "xy"            



set term x11
set output
#set term wxt
replot 
pause -1


set terminal postscript eps enhanced color
set output 'spec2.eps'
set title "down spin on (4,1)"

plot 'out.spec2.dat' using 1:17  w l title "3z^2-r^2",\
     'out.spec2.dat' using 1:18  w l title "zx",\
     'out.spec2.dat' using 1:19  w l title "yz",\
     'out.spec2.dat' using 1:20 w l title "x^2-y^2",\
     'out.spec2.dat' using 1:21 w l title "xy"            



set term x11
set output
#set term wxt
replot 
pause -1
