set xlabel "matrix order"
set ylabel "time (s)"
set title
#set ylabel "time (s)"
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
bb_bit (n) = 30.0;
bl_bit (n) = floor (log(67108864 / sqrt(1.0 * n)) / log (2));

bb_time (n) = 0.000001828748357*n*n*log(n) / bb_bit(n);
bl_time (n) = 0.0000000003732102412*n*n*n*log (n) / bl_bit (n);
lu_time (n) = 0.0000000003838888362*n*n*n / bl_bit (n);

#time per bit
bit_bb (t, n) = t / bb_bit(n);
bit_bl (t, n) = t / bl_bit (n);
bit_lu (t, n) = t / bl_bit (n);

plot [] []  "steps.table" using 2:(bit_bb($4, $2)) title "BB per bit" with points pt 1,\
			"steps.table" using 2:(bb_time($2)) title "Estimated BB" with lines lt 0,\
			"steps.table" using 2:(bit_bl($5, $2)) title "MINBLAS per bit" with points pt 3,\
			"steps.table" using 2:(bl_time($2)) title "Estimated BL" with lines lt 3,\
			"steps.table" using 2:(bit_lu($6, $2)) title "LUBLAS per bit" with points pt 6,\
			"steps.table" using 2:(lu_time($2)) title "Estimated LU" with lines lt -1 lw 2

#plot [] []  "time_table" using 2:4 title "BB" with points pt 1,\
#			"time_table" using 2:(bb_time($2)) title "Predicted BB" with linespoint pt 2,\
#			"time_table" using 2:6 title "MINBLAS" with points pt 3,\
#			"time_table" using 2:(bl_time($2)) title "Predicted BL" with linespoint pt 4,\
#			"time_table" using 2:8 title "LUBLAS" with points pt 6,\
#			"time_table" using 2:(lu_time($2)) title "Predicted LU" with linespoint pt 8

#set terminal postscript eps enhanced color solid
set terminal postscript eps enhanced monochrome solid "Helvetica" 20
set output "steps.eps"
replot
set terminal x11
replot
