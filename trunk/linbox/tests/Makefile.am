SUBDIRS = data

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) 
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

TESTS = 					\
	test-vector-domain test-sparse-matrix	\
	test-diagonal test-minpoly		\
	test-inverse test-large-modular

noinst_PROGRAMS = $(TESTS)

test_vector_domain_SOURCES = 		\
	test-vector-domain.C		\
	test-common.C	test-common.h

test_sparse_matrix_SOURCES = 		\
	test-sparse-matrix.C		\
	test-common.C	test-common.h

test_diagonal_SOURCES = 		\
	test-diagonal.C			\
	test-common.C	test-common.h

test_minpoly_SOURCES = 			\
	test-minpoly.C			\
	test-common.C	test-common.h

test_inverse_SOURCES =			\
	test-inverse.C			\
	test-common.C	test-common.h

test_large_modular_SOURCES =			    \
	test-large-modular.C	test-field-common.h \
	test-common.C	test-common.h

test_sparse_matrix_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la
test_vector_domain_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la
test_diagonal_DEPENDENCIES      = $(top_builddir)/linbox/liblinbox.la
test_minpoly_DEPENDENCIES       = $(top_builddir)/linbox/liblinbox.la
test_inverse_DEPENDENCIES       = $(top_builddir)/linbox/liblinbox.la
test_large_modular_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la
#	test-rational-solver-adaptive.C \
#	test-common.h

#test_smith_form_adaptive_SOURCES=    \
#	test-common.h	\
#	test-smith-form-adaptive.C
