INCLUDES=-I$(top_srcdir)
#LDFLAGS= -z muldefs
AM_CXXFLAGS = -g -DLinBoxSrcOnly -Wall 
AM_CFLAGS = -g -Wall

SUBDIRS = data

AM_CPPFLAGS=-g $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(ATLAS_CFLAGS) 
LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(ATLAS_LIBS) 



# Put an entry for your test-mycomponent.C here.  
# Don't forget to add a test_mycomponent_SOURCES entry below.
BASIC_TESTS =  			\
	test-commentator	\
	test-gmp-rational	\
	test-unparametric-field	\
	test-randiter-nonzero	\
	test-subvector		\
	test-vector-domain	\
	test-matrix-domain	\
	test-dense		\
	test-submatrix		\
	test-subiterator	\
	test-sparse		\
	test-scalar-matrix	\
	test-butterfly		\
	test-minpoly		\
	test-diagonal		\
	test-dif		\
	test-trace		\
	test-rank		\
	test-inverse		\
	test-moore-penrose	\
	test-param-fuzzy        \
	test-hom	        \
	test-image-field        \
	test-triplesbb		\
	test-direct-sum		\
	test-companion		\
	test-hilbert		\
	test-modular		\
	test-modular-int	\
	test-modular-short	\
	test-modular-double	\
	test-modular-byte	\
	test-modular-balance-int\
	test-frobenius          \
	test-det		\
	test-zero-one 		\
	test-lu			\
	test-bitonic-sort


# doesn't compile
#	test-lu			\
#	test-solve		\
#	test-gf2		\
#
# compile with warnings but test fails
#	test-gssv_rank		\
#
# a benchmarker, not to be included in check.
#	test-fields		\
#
# historical tests, no longer in use
#	test-nag-sparse		\ replaced by tripplesbb


# no explicit test needed, I guess...
#	Transpose is tested in test-triplesbb
#	Compose is tested in test-butterfly

# removed test
#	test-nag-sparse	

GIVARO_TESTS =			\
	test-givaro-zpz	


NTL_TESTS = 			\
	test-ntl-zz_p		\
	test-ntl-ZZ_p		\
	test-ntl-RR		\
	test-ntl-toeplitz	\
	test-ntl-hankel   	\
	test-ntl-sylvester	\
	test-sum		\
	test-gssv_rank		\
	test-local-smith

LIDIA_TESTS =			\
	test-lidia-gfq

ATLAS_TESTS = 			\
	test-blas-domain	\
	test-fflapack		\
	test-iliopoulos-elimination \
	test-smith-form         \
	test-last-invariant-factor\
	test-rational-solver

if LINBOX_HAVE_NTL
USE_NTL_TESTS = $(NTL_TESTS)
endif

if LINBOX_HAVE_GIVARO
USE_GIVARO_TESTS = $(GIVARO_TESTS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_TESTS = $(LIDIA_TESTS)
endif

if LINBOX_HAVE_ATLAS
USE_ATLAS_TESTS = $(ATLAS_TESTS)
endif

TESTS =				\
	$(BASIC_TESTS)		\
	$(USE_GIVARO_TESTS)	\
	$(USE_NTL_TESTS)	\
	$(USE_LIDIA_TESTS)	\
	$(USE_ATLAS_TESTS)

CLEANFILES = $(TESTS)

EXTRA_PROGRAMS = \
	$(BASIC_TESTS)		\
	$(GIVARO_TESTS)		\
	$(NTL_TESTS)		\
	$(ATLAS_TESTS)		\
	$(LIDIA_TESTS)		


test_commentator_SOURCES =		\
	test-commentator.C		\
	test-common.h

test_modular_SOURCES =			\
	test-modular.C			\
	test-generic.h			\
	test-common.h

test_hom_SOURCES =			\
	test-hom.C			\
	test-generic.h			\
	test-common.h

test_image_field_SOURCES =		\
	test-image-field.C		\
	test-generic.h			\
	test-common.h

#test_gf2_SOURCES =			\
#	test-gf2.C			\
#	test-generic.h			\
#	test-common.h
#
test_gmp_rational_SOURCES =		\
	test-gmp-rational.C		\
	test-generic.h			\
	test-common.h

test_unparametric_field_SOURCES =	\
	test-unparametric-field.C	\
	test-generic.h			\
	test-common.h

test_randiter_nonzero_SOURCES =		\
	test-randiter-nonzero.C		\
	test-common.h

test_subiterator_SOURCES = 		\
	test-subiterator.C		\
	test-common.h

test_subvector_SOURCES = 		\
	test-subvector.C		\
	test-common.h

test_vector_domain_SOURCES = 		\
	test-vector-domain.C		\
	test-common.h

test_matrix_domain_SOURCES = 		\
	test-matrix-domain.C		\
	test-common.h

test_dense_SOURCES = 			\
	test-dense.C			\
	test-common.h

test_sparse_SOURCES = 			\
	test-sparse.C			\
	test-common.h

test_submatrix_SOURCES = 		\
	test-submatrix.C		\
	test-common.h

test_scalar_matrix_SOURCES = 		\
	test-scalar-matrix.C		\
	test-generic.h			\
	test-common.h

test_butterfly_SOURCES = 		\
	test-butterfly.C		\
	test-generic.h			\
	test-common.h

test_minpoly_SOURCES = 			\
	test-minpoly.C			\
	test-common.h

test_diagonal_SOURCES = 		\
	test-diagonal.C			\
	test-generic.h			\
	test-common.h

test_sum_SOURCES = 			\
	test-sum.C			\
	test-generic.h			\
	test-common.h

test_dif_SOURCES = 			\
	test-dif.C	test-generic.h	\
	test-common.h

test_triplesbb_SOURCES =		\
	test-triplesbb.C		\
	test-generic.h			\
	test-common.h

test_hilbert_SOURCES =			\
	test-hilbert.C	test-generic.h	\
	test-common.h

test_ntl_toeplitz_SOURCES =		\
	test-ntl-toeplitz.C		\
	test-common.h

test_ntl_hankel_SOURCES =		\
	test-ntl-hankel.C		\
	test-common.h


test_ntl_sylvester_SOURCES =		\
	test-ntl-sylvester.C		\
	test-common.h

test_trace_SOURCES = 			\
	test-trace.C			\
	test-common.h

test_det_SOURCES = 			\
	test-det.C			\
	test-common.h

test_rank_SOURCES =			\
	test-rank.C			\
	test-common.h

#test_solve_SOURCES =			\
#	test-solve.C			\
#	test-common.h


test_direct_sum_SOURCES =		\
	test-direct-sum.C		\
	test-generic.h			\
	test-common.h

test_frobenius_SOURCES =		\
	test-frobenius.C		\
	test-generic.h			\
	test-common.h

test_companion_SOURCES =		\
	test-companion.C		\
	test-generic.h			\
	test-common.h

test_inverse_SOURCES =			\
	test-inverse.C			\
	test-common.h

test_moore_penrose_SOURCES =		\
	test-moore-penrose.C		\
	test-common.h

test_ntl_RR_SOURCES =			\
	test-ntl-RR.cpp			\
	test-generic.h			\
	test-common.h

test_ntl_ZZ_p_SOURCES =			\
	test-ntl-ZZ_p.cpp		\
	test-generic.h			\
	test-common.h

test_ntl_zz_p_SOURCES =			\
	test-ntl-zz_p.cpp		\
	test-generic.h			\
	test-common.h

test_givaro_zpz_SOURCES =		\
	test-givaro-zpz.C		\
	test-generic.h			\
	test-common.h

test_lidia_gfq_SOURCES =		\
	test-lidia-gfq.C		\
	test-generic.h			\
	test-common.h

test_gssv_rank_SOURCES =	       \
	test-gssv_rank.C	       

test_param_fuzzy_SOURCES =	       \
	test-param-fuzzy.C	       \
	test-generic.h		       \
	test-common.h

test_lu_SOURCES =                      \
	test-lu.C			 

test_zero_one_SOURCES =			\
	test-zero-one.C			\
	test-generic.h			\
	test-common.h

test_modular_int_SOURCES =		\
	test-modular-int.C 		\
	test-common.h	

test_modular_short_SOURCES =		\
	test-modular-short.C		\
	test-generic.h			\
	test-common.h

test_modular_double_SOURCES =		\
	test-modular-double.C		\
	test-generic.h			\
	test-common.h

test_modular_byte_SOURCES =		\
	test-modular-byte.C		\
	test-generic.h			\
	test-common.h

test_modular_balance_int_SOURCES =             \
        test-modular-balance-int.C             \
        test-generic.h                  \
        test-common.h

test_fields_SOURCES =			\
	test-fields.C

test_rational_solver_SOURCES =		\
	test-common.h	\
	test-rational-solver.C

test_iliopoulos_elimination_SOURCES=    \
	test-common.h	\
	test-iliopoulos-elimination.C

test_smith_form_SOURCES=    \
	test-common.h	\
	test-smith-form.C

test_last_invariant_factor_SOURCES=    \
	test-common.h	\
	test-last-invariant-factor.C

test_local_smith_SOURCES=    \
	test-common.h	\
	test-local-smith.C

test_bitonic_sort_SOURCES=    \
        test-common.h   \
        test-bitonic-sort.C

#test_sparse_rational_solver_SOURCES =		\
#	test-common.h	\
#	test-sparse-rational-solver.C

test_blas_domain_SOURCES=		\
	test-common.h   \
	test-blas-domain.C

test_fflapack_SOURCES=		\
	test-common.h   \
	test-fflapack.C
#	test-rational-solver-adaptive.C \
#	test-common.h

#test_smith_form_adaptive_SOURCES=    \
#	test-common.h	\
#	test-smith-form-adaptive.C
