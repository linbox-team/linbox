INCLUDES=-I$(top_srcdir)

SUBDIRS = data

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
LDADD = $(top_builddir)/linbox/liblinbox.la $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS)

# Put an entry for your test-mycomponent.C here.  
# Don't forget to add a test_mycomponent_SOURCES entry below.
BASIC_TESTS = 			\
	test-modular		\
	test-unparametric-field	\
	test-randiter-nonzero	\
	test-subvector		\
	test-vector-domain	\
	test-dense-matrix	\
	test-dense		\
	test-submatrix		\
	test-subiterator	\
	test-sparse		\
	test-scalar-matrix	\
	test-minpoly		\
	test-diagonal		\
	test-sum		\
	test-trace		\
	test-det		\
	test-rank		\
	test-solve		\
	test-inverse		\
	test-moore-penrose	\
	test-param-fuzzy

# known NTL_TESTS (1Aug02) 
#../macros/ntl-check.m4:NTL_TESTS="test-ntl-zz_p test-ntl-ZZ_p test-ntl-RR"

# known GIVARO_TESTS (1Aug02) 
#../macros/givaro-check.m4:GIVARO_TESTS="test-givaro-zpz"

# known LIDIA_TESTS (1Aug02) 
#../macros/lidia-check.m4:LIDIA_TESTS="test-lidia-gfq"

TESTS =				\
	$(BASIC_TESTS)		\
	$(NTL_TESTS)		\
	$(GIVARO_TESTS)         \
	$(LIDIA_TESTS)		\
	$(SOLVE_TESTS)

EXTRA_PROGRAMS = $(BASIC_TESTS) test-ntl-RR test-ntl-ZZ_p test-ntl-zz_p test-givaro-zpz test-lidia-gfq test-toeplitz test-lu

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

test_subiterator_SOURCES = 		\
	test-subiterator.C		\
	test-common.C	test-common.h

test_subvector_SOURCES = 		\
	test-subvector.C		\
	test-common.C	test-common.h

test_vector_domain_SOURCES = 		\
	test-vector-domain.C		\
	test-common.C	test-common.h

test_dense_matrix_SOURCES = 		\
	test-dense-matrix.C		\
	test-common.C	test-common.h

test_dense_SOURCES = 			\
	test-dense.C			\
	test-common.C	test-common.h

test_sparse_SOURCES = 			\
	test-sparse.C			\
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

test_toeplitz_SOURCES =			\
	test-toeplitz.C			\
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

test_solve_SOURCES =			\
	test-solve.C			\
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

test_givaro_zpz_SOURCES =		\
	test-givaro-zpz.C		\
	test-generic.h			\
	test-common.C	test-common.h

test_lidia_gfq_SOURCES =		\
	test-lidia-gfq.C		\
	test-generic.h			\
	test-common.C	test-common.h

test_lu_SOURCES =                      \
        test-lu.h                      \
	test-lu.C

test_param_fuzzy_SOURCES =			\
	test-param-fuzzy.C			\
	test-generic.h			\
	test-common.C	test-common.h

#test_stub_SOURCES =			\
#	test-stub.C			\
#	test-generic.h			\
#	test-common.C	test-common.h
#	test-rational-solver-adaptive.C \
#	test-common.h

#test_smith_form_adaptive_SOURCES=    \
#	test-common.h	\
#	test-smith-form-adaptive.C
