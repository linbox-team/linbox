TESTS=rational_rank_test
#
#noinst_PROGRAMS=rank_sparse_reord_gauss
#
#rank_sparse_reord_gauss_SOURCES=rank_sparse_reord_gauss.C
#
#INCLUDES=-I$(prefix)/include

rational_rank_test: rank_sparse_reord_gauss
	@echo "./rank_sparse_reord_gauss s=3mat.gz" >rational_rank_test
	@chmod +x rational_rank_test

%:%.C
	$(CXX) $(OPTFLAGS) $< -o $@ -I$(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS)
