set xlabel "matrix order"
set ylabel "time (s)"
set title
#set title  "relative efficiency"
set key bottom right
set logscale y 
set logscale x
set pointsize 2
#show pointsize

#plot [] [] "time_table" using 1:($6 / $2) title "MINBLAS" with points pointtype 3 pointsize 2, "time_table" using 1:($8 / $2) title "LUBLAS" with linespoint pt 1 ps 2

#plot [] [] "time_table" using 2:($4 / $6) title "MINBLAS" with linespoint pt 1 ps 2, \
#"time_table" using 2:($4 / $8) title "LUBLAS" with linespoint pt 1 ps 2

#formula (x, y) = y;
bb_time (n) = 0.002240174633*n*n  + 0.00006676817979*n*n*n;
bl_time (n) = 0.00001285796899*n*n*n  + 0.00000002885270500*n*n*n*n;
lu_time (n) = 0.000009645194584*n*n*n  + 0.0000000002217193167*n*n*n*n;

plot [][] \
			"predict.table" using 2:(bb_time($2)) title "Predicted BB" with lines lt 0,\
			"predict.table" using 2:(bl_time($2)) title "Predicted BL" with lines lt 1,\
			"predict.table" using 2:(lu_time($2)) title "Predicted LU" with lines lt -1 lw 2

#set terminal postscript eps enhanced color solid
set terminal postscript eps enhanced monochrome solid "Helvetica" 20
set output "predict.eps"
replot
set terminal x11
replot
