SUBDIRS=Gauss integertest


## FIXME this should not be here, but I don't know where it should be.
#env:
#	set -a LD_LIBRARY_PATH=${gmp_prefix}/lib:${LD_LIBRARY_PATH}
#	set -a LD_LIBRARY_PATH=${ntl_prefix}/lib:${LD_LIBRARY_PATH}
#	export LD_LIBRARY_PATH
#	echo ${LD_LIBRARY_PATH}
#

## When I compile test.cpp with the following command line
## it has type problems, but it seems to find all the header files at least.
test: test.cpp test.h test_base.h test_linbox.h
	g++ test.cpp -I $(prefix)/include $(NTL_CFLAGS) $(GMP_CFLAGS) $(GMP_LIBS) $(NTL_LIBS) -DNTL_STD_CXX -o test

