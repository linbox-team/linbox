LBFLAGS=-enable-shared -DLinBoxSrcOnly -DHAVE_CONFIG_H -D__LINBOX_HAVE_CBLAS -DDISABLE_COMMENTATOR -DBIG
LDPLACES=-lstdc++ -lm -lc -lgcc_s -lntl /usr/lib/libgmpxx.so /usr/lib/libgmp.so /usr/local/lib/libgivaro.so -L/usr/lib/sse2 -latlas -lblas -L/usr/lib/sse2/atlas -llapack

#-lstdc++ -lm -lc -lgcc_s -lntl /usr/lib/libgmpxx.so /usr/lib/libgmp.so /usr/lib/libgivaro.so -L/usr/lib/sse2/ -latlas -L/usr/lib/sse2/atlas -llapack -lcblas 

liblinboxfunctions.so: linboxfunctions.o linboxfunctions_wrap.o
	gcc -shared linboxfunctions.o linboxfunctions_wrap.o -o liblinboxfunctions.so ${LDPLACES}

linboxfunctions.o: linboxfunctions.C 
	gcc -fpic ${LBFLAGS} -c linboxfunctions.C -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux  -I/home/fendt/linbox/

test-smith.o: test-smith.C linboxfunctions.C makefile 
	g++ -fpic ${LBFLAGS} -c test-smith.C -I/home/fendt/linbox/

test-smith: test-smith.o  makefile
	g++ -shared test-smith.o ${LDPLACES}

linboxfunctions_wrap.o: linboxfunctions_wrap.cxx 
	gcc -fpic ${LBFLAGS} -c linboxfunctions_wrap.cxx -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux  -I/home/fendt/linbox/

clean: 
	rm *.o *.so *.a
