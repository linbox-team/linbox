set xlabel "matrix order" 
set ylabel "bit length"
set title
#set title  "relative efficiency"
set key bottom right
show key
#set logscale y 
set logscale x
set pointsize 2
#show pointsize

#plot [] [] "time_table" using 1:($6 / $2) title "MINBLAS" with points pointtype 3 pointsize 2, "time_table" using 1:($8 / $2) title "LUBLAS" with linespoint pt 1 ps 2

#plot [] [] "time_table" using 2:($4 / $6) title "MINBLAS" with linespoint pt 1 ps 2, \
#"time_table" using 2:($4 / $8) title "LUBLAS" with linespoint pt 1 ps 2

#formula (x, y) = y;
bits (n) = 176.5334824 + 18.42786376 * log(n);

plot [] []  "bits.table" using 2:3 title "actual bit length" with points pt 1,\
			"bits.table" using 2:(bits($2)) title "estimated by: 176.5334824+18.42786376ln(n)" with lines

#set terminal postscript eps enhanced color solid
set terminal postscript eps enhanced monochrome solid "Helvetica" 20
set output "bits.eps"
replot
set terminal x11
replot
