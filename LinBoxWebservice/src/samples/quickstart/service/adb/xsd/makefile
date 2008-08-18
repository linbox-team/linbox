LBFLAGS=-enable-shared -DLinBoxSrcOnly -DHAVE_CONFIG_H -D__LINBOX_HAVE_CBLAS -DDISABLE_COMMENTATOR -DBIG
LDPLACES=-llapack -lcblas -latlas -lntl /usr/lib/libgmpxx.so /usr/lib/libgmp.so -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/../../.. -lstdc++ -lm -lc -lgcc_s /usr/local/lib/libgivaro.a

liblinboxfunctions.so: linboxfunctions.o linboxfunctions_wrap.o
	gcc -shared linboxfunctions.o linboxfunctions_wrap.o -o liblinboxfunctions.so ${LDPLACES}

linboxfunctions.o: linboxfunctions.C 
	gcc -fpic ${LBFLAGS} -c linboxfunctions.C -I/opt/jdk1.5.0_05/include -I/opt/jdk1.5.0_05/include/linux  -I/home/fendt/linbox/

linboxfunctions_wrap.o: linboxfunctions_wrap.cxx 
	gcc -fpic ${LBFLAGS} -c linboxfunctions_wrap.cxx -I/opt/jdk1.5.0_05/include -I/opt/jdk1.5.0_05/include/linux  -I/home/fendt/linbox/

