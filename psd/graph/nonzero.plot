set xlabel "matrix order"
set ylabel "number of non-zero entries"
#set title  "relative efficiency"
set title
set key bottom right
set logscale y 
set logscale x
set pointsize 2
#show pointsize

#plot [] [] "time_table" using 1:($6 / $2) title "MINBLAS" with points pointtype 3 pointsize 2, "time_table" using 1:($8 / $2) title "LUBLAS" with linespoint pt 1 ps 2

#plot [] [] "time_table" using 2:($4 / $6) title "MINBLAS" with linespoint pt 1 ps 2, \
#"time_table" using 2:($4 / $8) title "LUBLAS" with linespoint pt 1 ps 2

#formula (x, y) = x * x * x * x * log(x) / y / 1000000;
#formula (x, y) = y;



formula (n) = 87.7 * n + 20.0 * n * log (n);

#plot [] []  "nonzero_table" using 2:3 title "NonZeroEntries" with linespoint 1,\
#"nonzero_table" using 2:(formula($2)) title "Estimated" with linespoint 2

#plot [] []  "nonzero_table" using 2:(formula2($2, $3)) title "relative errors: 87.7n + 20nlogn" with linespoint 1

plot [] []  "nonzero.table" using 2:(formula($2)) title "estimated by: 87.7n + 20nlogn" with lines 1,\
			"nonzero.table" using 2:3 title "Actual data" with points pt 2
			
set terminal postscript eps enhanced monochrome solid "Helvetica" 20
set output "nonzero.eps"
replot
set terminal x11
replot
