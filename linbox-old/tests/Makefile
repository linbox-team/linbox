# this is the oldstyle Makefile so will likely disappear.
# Save this info somewhere:

## When I compile test.cpp with the following command line
## it has type problems, but it seems to find all the header files at least.
#g++ test.cpp -I /usa/saunders/new/linbox/include/ -I /usa/saunders/algebra/NTL/include/ -DNTL_STD_CXX
# end of -bds-jul01 preamble 

# integertry works with integer.h, not yet with lin_int_giv.h (because of gmp.h)
LBinclude=/usa/saunders/new/linbox/include/
GMPinclude=/usr/local/algebra/gmp/include/
integertry: integertry.cc
	g++ integertry.cc -I $(LBinclude) -I $(GMPinclude)

TESTDIRS=Random Gauss Sequence
#include Makefile.rules

tests:
	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)

clean:
	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)

distclean:
	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)


