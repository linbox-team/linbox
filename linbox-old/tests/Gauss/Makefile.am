TESTS=rank_sparse_reord_gauss
#
#noinst_PROGRAMS=rank_sparse_reord_gauss
#
#rank_sparse_reord_gauss_SOURCES=rank_sparse_reord_gauss.C
#
#INCLUDES=-I$(prefix)/include

%:%.C
	$(CXX) $< -o $@ -I$(prefix)/include $(GMP_CFLAGS) $(GMP_LIBS)
