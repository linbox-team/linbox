# To make some test programs

TESTS=integertry

## When I compile test.cpp with the following command line
## it has type problems, but it seems to find all the header files at least.
test: test.cpp test.h test_base.h test_linbox.h
	g++ test.cpp -I $(prefix)/include $(NTL_CFLAGS) $(NTL_LIBS) -DNTL_STD_CXX -o test

# integertry works with lin_integer.h 
integertry: integertry.cc
	g++ integertry.cc -I $(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS) -o integertry

SUBDIRS=Gauss
