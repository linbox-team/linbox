SUBDIRS=

TESTS= dotprod vaxpy sum scalarprod

%:%.cpp
	$(CXX) $(OPTFLAGS) $< -o $@ -O -I$(prefix)/include \
	$(GMP_CFLAGS) $(GMP_LIBS) \
	$(NTL_CFLAGS) $(NTL_LIBS) 

#
#integercheck: integercheck.cc
#	g++ integercheck.cc -I/usr/local/algebra/gmp-4.0/include -L/usr/local/algebra/gmp-4.0/lib -R/usr/local/algebra/gmp-4.0/lib -R/usr/local/gnu/gcc/lib -lgmp
#
#	
clean: 
	rm $(TESTS)
