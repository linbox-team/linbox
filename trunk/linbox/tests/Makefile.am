SUBDIRS = data

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) 
LDADD = $(GMP_LIBS) $(NTL_LIBS) $(top_builddir)/linbox/liblinbox.la

TESTS = test-sparse-matrix
noinst_PROGRAMS = test-sparse-matrix

test_sparse_matrix_SOURCES = \
	test-sparse-matrix.C		\
	test-common.C	test-common.h

test_sparse_matrix_DEPENDENCIES=$(top_builddir)/linbox/liblinbox.la
#	test-rational-solver-adaptive.C \
#	test-common.h

#test_smith_form_adaptive_SOURCES=    \
#	test-common.h	\
#	test-smith-form-adaptive.C
