SUBDIRS=Gauss integertest 

TESTS=integertry test toeptest test_snf

#vpath %.cpp $(shell find . -type d)
#
#
#Field_Padic double float modular fuzzy param_modular param_fuzzy ntl_ZZ_p ntl_zz_p ntl_rr double_envelope float_envelope modular_envelope fuzzy_envelope param_modular_envelope param_fuzzy_envelope ntl_ZZ_p_envelope ntl_zz_p_envelope ntl_rr_envelope abstract_double abstract_float abstract_modular abstract_fuzzy abstract_param_modular abstract_param_fuzzy test:	%:	%.cpp;
#	g++ $(CXXFLAGS) -I$(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS) $(NTL_CFLAGS) $(NTL_LIBS) -o $@ $(filter %.cpp %.o, $^) -lgmp -lntl
#

%:%.cc
	$(CXX) $(OPTFLAGS) $< -o $@ -I$(prefix)/include \
	$(GMP_CFLAGS) $(GMP_LIBS) \
	$(NTL_CFLAGS) $(NTL_LIBS) 

