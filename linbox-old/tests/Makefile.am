# To make some test programs

## When I compile test.cpp with the following command line
## it has type problems, but it seems to find all the header files at least.
test: test.cpp test.h test_base.h test_linbox.h
	g++ test.cpp -I $(prefix)/include $(NTL_CFLAGS) $(NTL_LIBS) -DNTL_STD_CXX

# integertry works with integer.h, not yet with lin_integer.h (io prob)
integertry: integertry.cc
	g++ integertry.cc -I $(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS)

#SUBDIRS=

## Save this oldstyle Makefile info as template.  
## We might follow this info in running the subdir tests
#
#TESTDIRS=Random Gauss Sequence
##include Makefile.rules
#
#tests:
#	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)
#
#clean:
#	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)
#
#distclean:
#	$(foreach SDIR,$(TESTDIRS),$(MAKE) -C $(SDIR) $@;)
