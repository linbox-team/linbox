
EXTRA_DIST= ex_field.C ex_field_archetype.C ex_field_wrapper.C

%:%.C
	$(CXX) $(OPTFLAGS) $< -o $@ -I$(prefix)/include $(GMP_CFLAGS) $(GIVARO_CFLAGS) $(NTL_CFLAGS) $(GMP_LIBS) $(GIVARO_LIBS) $(NTL_LIBS) -lgmp -lgivaro -lntl
