SUBDIRS=

TESTS=	wiedemann_det wiedemann_minpoly wiedemann_linsolve1

%:%.cpp
	$(CXX) $(OPTFLAGS) $< -o $@ -O -I$(prefix)/include \
	$(GMP_CFLAGS) $(GMP_LIBS) \
	$(NTL_CFLAGS) $(NTL_LIBS) 

clean: 
	rm $(TESTS)
