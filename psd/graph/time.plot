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

nano = 1000000000;

bb_time (n) = 10049.57636*n*n*n*log (n) / nano;
bl_time (n) = 6.075421680*n*n*n*n*log(n) / nano;
lu_time (n) = 10.0234965*n*n*n*n / nano;

set linestyle 1 lt 1 lw 1
set linestyle 2 lt 2 lw 1 
set linestyle 3 lt 3 lw 1 
set linestyle 4 lt 6 lw 1 
set linestyle 5 lt 1 lw 3 
set linestyle 6 lt 2 lw 3 
set linestyle 7 lt 3 lw 3 
set linestyle 8 lt 6 lw 3

plot [] []  "time.table" using 2:4 title "BB time" with points pt 1,\
			"time.table" using 2:(bb_time($2)) title "Estimated BB" with lines lt 0,\
			"time.table" using 2:6 title "MINBLAS timeS" with points pt 3,\
			"time.table" using 2:(bl_time($2)) title "Estimated MINBLAS" with lines lt 1,\
			"time.table" using 2:8 title "LUBLAS time" with points pt 6,\
			"time.table" using 2:(lu_time($2)) title "Estimated LUBLAS" with lines lt -1 lw 2

#plot [] []  "time_table" using 2:4 title "BB" with points pt 1,\
#			"time_table" using 2:(bb_time($2)) title "Predicted BB" with linespoint pt 2,\
#			"time_table" using 2:6 title "MINBLAS" with points pt 3,\
#			"time_table" using 2:(bl_time($2)) title "Predicted BL" with linespoint pt 4,\
#			"time_table" using 2:8 title "LUBLAS" with points pt 6,\
#			"time_table" using 2:(lu_time($2)) title "Predicted LU" with linespoint pt 8

#set terminal postscript eps enhanced color solid
set terminal postscript eps enhanced monochrome solid "Helvetica" 20
set output "time.eps"
replot
set terminal x11
replot
