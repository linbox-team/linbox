SUBDIRS=

TESTS= diagonal hilbert sparsemat compose butterfly transpose

%:%.cpp
	$(CXX) $(OPTFLAGS) $< -o $@ -O -I$(prefix)/include \
	$(GMP_CFLAGS) $(GMP_LIBS) \
	$(NTL_CFLAGS) $(NTL_LIBS) 

clean: 
	rm $(TESTS)
