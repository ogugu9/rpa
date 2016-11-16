set terminal postscript eps enhanced color
set output 'chi0.eps'

set size square

#kpi is the # of k in pi
k2pi = 120
kpi = k2pi / 2
kpi = 1

max=200  # 強度の最大値
min=0	 # 強度の最小値
#set palette defined (min "black", min "blue", 0 "white", max/2.0 "orange", max "red")

#set palette defined (0 "blue",  max/2.0 "green", max "red")
#set palette defined (-40 "violet", -5 "blue", 0 "dark-green", 5 "yellow", 10 "red")
set palette defined (0 "dark-green", 20 "yellow", 100 "red")
set pm3d corners2color c1



#set cbrange[min:max]  # 強度の範囲

#set xrange [0:40]
#set yrange[0:40]
#set zrange[min:max]
set pm3d map
#set pm3d at bs
set tics out
splot 'out.chi.dat' using 1:2:3



set terminal x11
set output
replot
pause -1

