set xlabel "matrix order"
set ylabel "bits"
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

bb_b (n) = 30;
blas_b (n) = floor (log (67108864 / log (n)) / log (2.0));

bb_bit (n, p) = p * bb_b(n);
blas_bit (n, p) = p * blas_b (n);

estimated_bit (n) = 150.1 * n + 17.30 * n * log(n)

plot [] []  "bit_required.table" using 2:(bb_bit($2, $5)) title "Actual required" with points pt 1,\
			"bit_required.table" using 2:(estimated_bit($2)) title "Estimated by:150.1n + 17.3nln(n) " with lines

#plot [] []  "time_table" using 2:4 title "BB" with points pt 1,\
#			"time_table" using 2:(bb_time($2)) title "Predicted BB" with linespoint pt 2,\
#			"time_table" using 2:6 title "MINBLAS" with points pt 3,\
#			"time_table" using 2:(bl_time($2)) title "Predicted BL" with linespoint pt 4,\
#			"time_table" using 2:8 title "LUBLAS" with points pt 6,\
#			"time_table" using 2:(lu_time($2)) title "Predicted LU" with linespoint pt 8

#set terminal postscript eps enhanced color solid
set terminal postscript eps enhanced monochrome solid "Helvetica" 20
set output "bit_required.eps"
replot
set terminal x11
replot
