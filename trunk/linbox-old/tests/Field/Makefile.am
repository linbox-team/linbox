SUBDIRS=

TESTS=double envelope modular fuzzy param_modular param_fuzzy ntl_zz_p \
	abstract_double abstract_modular abstract_fuzzy \
	abstract_param_modular abstract_param_fuzzy

%:%.cpp	field.h
	$(CXX) $(OPTFLAGS) $< -o $@ -O -I$(prefix)/include \
	$(GMP_CFLAGS) $(GMP_LIBS) \
	$(NTL_CFLAGS) $(NTL_LIBS) 

clean: 
	rm $(TESTS)
