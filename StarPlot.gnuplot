#This gnuplot file should be called by the StarPolt.py
set multiplot 
set size 0.4,1
set origin 0,0
set pm3d
unset surface
set key outside
splot '/tmp/gridtoplot.table' notitle
unset surface
set view map
set contour
set key outside
set size 0.6,0.6
set origin 0.35,0.2
set size ratio 1
splot '/tmp/gridtoplot.table' notitle
unset multiplot
pause -1