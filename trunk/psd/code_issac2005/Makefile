all : psd-minpoly-bb psd-minpoly-lu psd-minpoly-blas

%:%.C
	g++ -O3 -DLinBoxSrcOnly -I /home/wan/linbox $< -lgmp -lgmpxx /home/wan/linbox/linbox/.libs/liblinbox.a -o $@ -lcblas -lclapack -latlas -lntl

nodense: stats.C
	g++ -O3 -DLinBoxSrcOnly -DNodense -I /home/wan/linbox stats.C -lgmp -lgmpxx /home/wan/linbox/linbox/.libs/liblinbox.a -o statsnodense -lcblas -lclapack -latlas -lntl
clean :
	rm -f psd-minpoly-bb psd-minpoly-lu psd-minpoly-blas

run : blas_out lu_out blackbox_out

blas_out : psd-minpoly-blas
	blas_run

lu_out : psd-minpoly-lu
	lu_run

blackbox_out : psd-minpoly-bb
	blackbox_run
