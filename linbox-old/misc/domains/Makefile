# Makefile for running the tests during 3/00  -bds

CPP = g++
# A C++ compiler, e.g., g++, supporting templates well

NTL = /usr/local/algebra/NTL
LIB = $(NTL)/src/ntl.a

CPPFLAGS =  -I$(NTL)/include -O2 -mv8

# Useful compile-time options (CFLAGS):
# -I../include -O2 -g -mv8 -qarch=ppc

a.out: svterm.cc
	$(CPP) $(CPPFLAGS) svterm.cc $(NTL)/src/ntl.a

.cc.o:
	$(CPP) $(CPPFLAGS) -c $< $(NTL)/src/ntl.a
	
Integers.o: Integers.cc Integers.h
	$(CPP) $(CPPFLAGS) Integers.cc -c 
	
test: test.cc Integers.o *.h
	$(CPP) $(CPPFLAGS) test.cc Integers.o $(NTL)/src/ntl.a -o test

vectest: useLB_vector.cc my_free_module.h LB_vector.h Integers.o
	$(CPP) $(CPPFLAGS) useLB_vector.cc Integers.o $(NTL)/src/ntl.a 
