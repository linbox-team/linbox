SUBDIRS=Gauss #integertest

TESTS=integertry test testminpoly toeptest

RUN=-R /usr/local/gnu/gcc/lib

vpath %.cpp $(shell find . -type d)


double float modular fuzzy param_modular param_fuzzy ntl_ZZ_p ntl_zz_p ntl_rr double_envelope float_envelope modular_envelope fuzzy_envelope param_modular_envelope param_fuzzy_envelope ntl_ZZ_p_envelope ntl_zz_p_envelope ntl_rr_envelope abstract_double abstract_float abstract_modular abstract_fuzzy abstract_param_modular abstract_param_fuzzy test:	%:	%.cpp;
	g++ $(CXXFLAGS) -I$(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS) $(NTL_CFLAGS) $(NTL_LIBS) -o $@ $(filter %.cpp %.o, $^) -lgmp -lntl

# integertry works with lin_integer.h 
integertry: integertry.cc
	g++ integertry.cc -I $(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS) -o integertry

testminpoly: testminpoly.cc 
	g++ testminpoly.cc $(GMP_CFLAGS) $(GMP_LIBS) $(NTL_CFLAGS) $(NTL_LIBS) -I $(prefix)/include $(RUN) -o testminpoly

toeptest: toeptest.cc 
	g++ toeptest.cc -I $(prefix)/include $(NTL_CFLAGS) $(NTL_LIBS) $(GMP_CFLAGS) $(GMP_LIBS) -o toeptest

fftToeptest: fftToeptest.cc fftToeplitz.h
	g++ fftToeptest.cc -I $(prefix)/include $(NTL_CFLAGS) $(NTL_LIBS) $(GMP_CFLAGS) $(GMP_LIBS) -R /usr/local/gnu/gcc/lib/ -o fftToeptest

