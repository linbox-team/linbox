
EXTRA_DIST= ex_fields.C ex_fields_archetype.C ex_fields_wrapper.C

%:%.C
	$(CXX) $(OPTFLAGS) $< -o $@ -I$(prefix)/include $(GMP_CFLAGS) $(GIVARO_CFLAGS) $(NTL_CFLAGS) $(GMP_LIBS) $(GIVARO_LIBS) $(NTL_LIBS) -lgmp -lgivaro -lntl
