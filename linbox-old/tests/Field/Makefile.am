SUBDIRS=

TESTS=double envelope modular fuzzy param_modular param_fuzzy ntl_zz_p

%:%.cpp
	$(CXX) $(OPTFLAGS) $< -o $@ -O -I$(prefix)/include \
	$(GMP_CFLAGS) $(GMP_LIBS) \
	$(NTL_CFLAGS) $(NTL_LIBS) 

clean: 
	rm $(TESTS)
