INCLUDES=-I$(top_srcdir)

SUBDIRS = data

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
LDADD = $(top_builddir)/linbox/liblinbox.la $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS)

#clobber non-compiling tests
SOLVE_TESTS=
#SOLVE_TESTS=test-solveNonsingular

TESTS =				\
	test-modular		\
	test-unparametric-field	\
	test-randiter-nonzero	\
	test-subvector		\
	test-vector-domain	\
	test-dense-matrix	\
	test-submatrix		\
	test-sparse0		\
	test-scalar-matrix	\
	test-minpoly		\
	test-diagonal		\
	test-sum		\
	test-trace		\
	test-det		\
	test-rank		\
	test-inverse		\
	test-moore-penrose	\
	$(NTL_TESTS)		\
	$(GIVARO_TESTS)         \
	$(LIDIA_TESTS)		\
	$(SOLVE_TESTS)		\
	test-dense-vector-domain


noinst_PROGRAMS = $(TESTS)
EXTRA_PROGRAMS = test-ntl-RR test-ntl-ZZ_p test-ntl-zz_p test-givaro-zpz test-lidia-gfq

test_modular_SOURCES =			\
	test-modular.C			\
	test-generic.h			\
	test-common.C	test-common.h

test_unparametric_field_SOURCES =	\
	test-unparametric-field.C	\
	test-generic.h			\
	test-common.C	test-common.h

test_randiter_nonzero_SOURCES =		\
	test-randiter-nonzero.C		\
	test-common.C	test-common.h

#test_subiterator_SOURCES = 		\
#	test-subiterator.C		\
#	test-common.C	test-common.h

test_subvector_SOURCES = 		\
	test-subvector.C		\
	test-common.C	test-common.h

test_vector_domain_SOURCES = 		\
	test-vector-domain.C		\
	test-common.C	test-common.h

test_dense_matrix_SOURCES = 		\
	test-dense-matrix.C		\
	test-common.C	test-common.h

#test_dense_matrix1_SOURCES = 		\
#	test-dense-matrix1.C		\
#	test-common.C	test-common.h
#
test_sparse0_SOURCES = 			\
	test-sparse0.C			\
	test-common.C	test-common.h

test_submatrix_SOURCES = 		\
	test-submatrix.C		\
	test-common.C	test-common.h

test_scalar_matrix_SOURCES = 		\
	test-scalar-matrix.C		\
	test-generic.h			\
	test-common.C	test-common.h

test_minpoly_SOURCES = 			\
	test-minpoly.C			\
	test-common.C	test-common.h

test_diagonal_SOURCES = 		\
	test-diagonal.C			\
	test-generic.h			\
	test-common.C	test-common.h

test_sum_SOURCES = 			\
	test-sum.C			\
	test-generic.h			\
	test-common.C	test-common.h

test_trace_SOURCES = 			\
	test-trace.C			\
	test-common.C	test-common.h

test_det_SOURCES = 			\
	test-det.C			\
	test-common.C	test-common.h

test_rank_SOURCES =			\
	test-rank.C			\
	test-common.C	test-common.h

test_inverse_SOURCES =			\
	test-inverse.C			\
	test-common.C	test-common.h

test_moore_penrose_SOURCES =		\
	test-moore-penrose.C		\
	test-common.C	test-common.h

test_ntl_RR_SOURCES =			\
	test-ntl-RR.cpp			\
	test-generic.h			\
	test-common.C	test-common.h

test_ntl_ZZ_p_SOURCES =			\
	test-ntl-ZZ_p.cpp		\
	test-generic.h			\
	test-common.C	test-common.h

test_ntl_zz_p_SOURCES =			\
	test-ntl-zz_p.cpp		\
	test-generic.h			\
	test-common.C	test-common.h

test_dense_vector_domain_SOURCES =      \
	test-dense-vector-domain.C      \
	test-common.C	test-common.h                   

test_givaro_zpz_SOURCES =		\
	test-givaro-zpz.C		\
	test-generic.h			\
	test-common.C	test-common.h

test_lidia_gfq_SOURCES =		\
	test-lidia-gfq.C		\
	test-generic.h			\
	test-common.C	test-common.h

#test_solveNonsingular_SOURCES =		\
#	test-solveNonsingular.C		\
#	test-common.C	test-common.h

test_modular_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_unparametric_field_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_randiter_nonzero_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_submatrix_DEPENDENCIES             = $(top_builddir)/linbox/liblinbox.la
test_vector_domain_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
#test_subiterator_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_sparse0_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_dense_matrix_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
#test_dense_matrix1_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
test_scalar_matrix_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_minpoly_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_diagonal_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_sum_DEPENDENCIES			= $(top_builddir)/linbox/liblinbox.la
test_det_DEPENDENCIES			= $(top_builddir)/linbox/liblinbox.la
test_rank_DEPENDENCIES			= $(top_builddir)/linbox/liblinbox.la
test_inverse_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_moore_penrose_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_ntl_RR_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_ntl_ZZ_p_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_ntl_zz_p_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_givaro_zpz_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
test_dense_vector_domain_DEPENDENCIES   = $(top_builddir)/linbox/liblinbox.la
test_lidia_gfq_DEPENDENCIES		= $(top_builddir)/linbox/liblinbox.la
#test_solveNonsingular_DEPENDENCIES	= $(top_builddir)/linbox/liblinbox.la
#	test-rational-solver-adaptive.C \
#	test-common.h

#test_smith_form_adaptive_SOURCES=    \
#	test-common.h	\
#	test-smith-form-adaptive.C
