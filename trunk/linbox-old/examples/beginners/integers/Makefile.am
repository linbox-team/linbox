
EXTRA_DIST= ex_integers.C

%:%.C
	$(CXX) $(OPTFLAGS) $< -o $@ -I$(prefix)/include $(GMP_CFLAGS) $(GIVARO_CFLAGS) $(NTL_CFLAGS) $(GMP_LIBS) $(GIVARO_LIBS) $(NTL_LIBS) 
