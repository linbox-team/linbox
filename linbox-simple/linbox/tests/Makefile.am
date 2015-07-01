# FACING_REALITY is used anywhere?
FACING_REALITY = 0
INCLUDES=-I$(top_srcdir)
AM_LDFLAGS = -llapack -lcblas -latlas 
AM_CXXFLAGS = -g -DLinBoxSrcOnly -Wall
#AM_CFLAGS = -g -Wall

SUBDIRS = data

AM_CPPFLAGS= $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)  
LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS)  

EXTRA_DIST=\
test-common.C

# Put an entry for your test-mycomponent.C here (or in one of the other test groups).  
# Don't forget to add a test_mycomponent_SOURCES entry below.
#
# Once present, tests should remain present, being moved from one group to another.
# Note that among the available groups are NONCOMPILING_TESTS and FAILING_TESTS, which
# will be ignored by "make check" and can be temporary storage for problematic tests.
BASIC_TESTS =  			\
        test-bitonic-sort$(EXEEXT)      \
	  test-blas-domain$(EXEEXT)\
        test-butterfly$(EXEEXT)         \
        test-commentator$(EXEEXT)       \
        test-companion$(EXEEXT)         \
        test-dense$(EXEEXT)             \
        test-det$(EXEEXT)               \
        test-diagonal$(EXEEXT)          \
        test-dif$(EXEEXT)               \
        test-direct-sum$(EXEEXT)        \
	  test-ffpack$(EXEEXT)	\
        test-frobenius$(EXEEXT)         \
        test-getentry$(EXEEXT)          \
        test-gf2$(EXEEXT)               \
        test-gmp-rational$(EXEEXT)      \
        test-hilbert$(EXEEXT)           \
        test-hom$(EXEEXT)               \
        test-inverse$(EXEEXT)           \
        test-isposdef$(EXEEXT)          \
        test-ispossemidef$(EXEEXT)      \
        test-lu$(EXEEXT)                \
        test-matrix-domain$(EXEEXT)     \
        test-matrix-stream$(EXEEXT)     \
        test-mg-block-lanczos$(EXEEXT)  \
        test-minpoly$(EXEEXT)           \
        test-modular$(EXEEXT)           \
        test-modular-balance-int\$(EXEEXT) \
	  test-modular-byte$(EXEEXT)      \
        test-modular-double$(EXEEXT)    \
        test-modular-int$(EXEEXT)       \
        test-modular-short$(EXEEXT)     \
        test-moore-penrose$(EXEEXT)     \
        test-randiter-nonzero$(EXEEXT)  \
        test-rank$(EXEEXT)              \
        test-scalar-matrix$(EXEEXT)     \
        test-solve$(EXEEXT)             \
        test-sparse$(EXEEXT)            \
        test-subiterator$(EXEEXT)       \
        test-submatrix$(EXEEXT)         \
        test-subvector$(EXEEXT)         \
        test-trace$(EXEEXT)             \
        test-triplesbb$(EXEEXT)         \
        test-vector-domain$(EXEEXT)     \
        test-zo$(EXEEXT)                \
        test-side-by-side$(EXEEXT)      \
        test-over-under$(EXEEXT)        \
        test-quad-matrix$(EXEEXT)       \
        test-zero-one$(EXEEXT)
# tests of components not adhering to LinBox requirements.
NONCOMPLIANT_TESTS =	\
	test-unparametric-field$(EXEEXT)\
	test-ntl-RR$(EXEEXT)	\
	test-param-fuzzy$(EXEEXT)
#These need clapack
#test-rational-solver-adaptive$(EXEEXT)#test-smith-form-adaptive$(EXEEXT) 

# this allows _SOURCES entry to remain
# except for "make dist"
NONEXISTING_TESTS =
#NONEXISTING_TESTS =		\
#	test-solutions$(EXEEXT)#
#test_solutions_SOURCES=    \
#	test-common.h	\
#	test-solutions.C
#

NONCOMPILING_TESTS =  		

FAILING_TESTS =		\
	test-la-block-lanczos$(EXEEXT)
USE_NONCOMPILING_TESTS = 
USE_FAILING_TESTS = 
# comment out these two lines to avoid facing reality
#USE_NONCOMPILING_TESTS += $(NONCOMPILING_TESTS)
#USE_FAILING_TESTS += $(FAILING_TESTS)

# a benchmarker, not to be included in check.
BENCHMARK_TESTS = 		\
	test-fields$(EXEEXT)	

# historical tests, no longer in use
#	test-nag-sparse$(EXEEXT)	  replaced by tripplesbb
#       test-lsp$(EXEEXT)	lsp.h superceded by factorized-matrix.h
# removed test
#	test-nag-sparse$(EXEEXT)#	test-image-field$(EXEEXT)      replaced by hom 
#       test-lsp$(EXEEXT)# I put test-gssv_rank$(EXEEXT)in hmrg:~saunders/gen_superlu.  It is not linbox until and unless it is made to work
# test-gssv$(EXEEXT)is an ntl test


# no explicit test needed, I guess...
#	Transpose is tested in test-triplesbb$(EXEEXT)#	Compose is tested in test-butterfly$(EXEEXT)
GIVARO_TESTS =			\
	test-charpoly$(EXEEXT)	\
	test-givaro-zpz$(EXEEXT)	\
	test-givaro-zpzuns$(EXEEXT)
NTL_TESTS = 				\
	test-ntl-lzz_p$(EXEEXT)		\
	test-ntl-ZZ_p$(EXEEXT)		\
	test-ntl-toeplitz$(EXEEXT)	\
	test-ntl-hankel$(EXEEXT)  		\
	test-ntl-sylvester$(EXEEXT)	\
	test-sum$(EXEEXT)		\
	test-smith-form-local$(EXEEXT)	\
	test-smith-form$(EXEEXT)		\
	test-smith-form-iliopoulos$(EXEEXT)\
	test-smith-form-adaptive$(EXEEXT)\
	test-smith-form-binary$(EXEEXT)	\
	test-rational-solver$(EXEEXT)	\
	test-rational-solver-adaptive$(EXEEXT)\
	test-last-invariant-factor$(EXEEXT)\
	test-toeplitz-det$(EXEEXT)
LIDIA_TESTS =			\
	test-lidia-gfq$(EXEEXT)
ATLAS_TESTS = 

if LINBOX_HAVE_NTL
USE_NTL_TESTS = $(NTL_TESTS)
endif

if LINBOX_HAVE_GIVARO
USE_GIVARO_TESTS = $(GIVARO_TESTS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_TESTS = $(LIDIA_TESTS)
endif


USE_ATLAS_TESTS = $(ATLAS_TESTS)


TESTS =				\
	$(BASIC_TESTS)		\
	$(USE_GIVARO_TESTS)	\
	$(USE_NTL_TESTS)	\
	$(USE_LIDIA_TESTS)	\
	$(USE_NONCOMPILING_TESTS)	\
	$(USE_FAILING_TESTS)	\
	$(USE_ATLAS_TESTS)


CLEANFILES = $(TESTS)

EXTRA_PROGRAMS = \
	$(BASIC_TESTS)		\
	$(GIVARO_TESTS)		\
	$(NTL_TESTS)		\
	$(ATLAS_TESTS)		\
	$(LIDIA_TESTS)		\
	$(BENCHMARK_TESTS) 		\
	$(NONEXISTING_TESTS)	\
	$(NONCOMPILING_TESTS)	\
	$(NONCOMPLIANT_TESTS)	\
	$(FAILING_TESTS)


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

test_gf2_SOURCES =			\
	test-gf2.C			\
	test-generic.h			\
	test-common.h

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

test_charpoly_SOURCES = 		\
	test-charpoly.C			\
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

test_toeplitz_det_SOURCES =		\
	test-toeplitz-det.C

test_trace_SOURCES = 			\
	test-trace.C			\
	test-common.h

test_getentry_SOURCES = 			\
	test-getentry.C			\
	test-common.h

test_minpoly_SOURCES = 			\
	test-minpoly.C			\
	test-common.h

test_det_SOURCES = 			\
	test-det.C			\
	test-common.h

test_rank_SOURCES =			\
	test-rank.C			\
	test-common.h

test_solve_SOURCES =			\
	test-solve.C			\
	test-common.h


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
	test-ntl-RR.C			\
	test-generic.h			\
	test-common.h

test_ntl_ZZ_p_SOURCES =			\
	test-ntl-ZZ_p.C			\
	test-generic.h			\
	test-common.h

test_ntl_lzz_p_SOURCES =		\
	test-ntl-lzz_p.C		\
	test-generic.h			\
	test-common.h

test_givaro_zpz_SOURCES =		\
	test-givaro-zpz.C		\
	test-generic.h			\
	test-common.h

test_givaro_zpzuns_SOURCES =		\
	test-givaro-zpzuns.C		\
	test-generic.h			\
	test-common.h

test_lidia_gfq_SOURCES =		\
	test-lidia-gfq.C		\
	test-generic.h			\
	test-common.h

test_param_fuzzy_SOURCES =	        \
	test-param-fuzzy.C	        \
	test-generic.h		        \
	test-common.h

test_lu_SOURCES =                       \
	test-lu.C			 

test_side_by_side_SOURCES =		\
	test-side-by-side.C		\
	test-generic.h			\
	test-common.h

test_over_under_SOURCES =		\
	test-over-under.C		\
	test-generic.h			\
	test-common.h

test_quad_matrix_SOURCES =		\
	test-quad-matrix.C		\
	test-generic.h			\
	test-common.h

test_zo_SOURCES =			\
	test-zo.C			\
	test-generic.h			\
	test-common.h

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

test_smith_form_iliopoulos_SOURCES=    \
	test-common.h	\
	test-smith-form-iliopoulos.C

test_smith_form_binary_SOURCES=    \
	test-common.h	\
	test-smith-form-binary.C

test_smith_form_SOURCES=    \
	test-common.h	\
	test-smith-form.C

test_last_invariant_factor_SOURCES=    \
	test-common.h	\
	test-last-invariant-factor.C

test_smith_form_local_SOURCES=    \
	test-common.h	\
	test-smith-form-local.C

test_bitonic_sort_SOURCES=    \
        test-common.h   \
        test-bitonic-sort.C

test_rational_solver_SOURCES =		\
	test-common.h	\
	test-rational-solver.C

test_blas_domain_SOURCES=		\
	test-common.h   \
	test-blas-domain.C

test_ffpack_SOURCES=		\
	test-common.h   \
	test-ffpack.C

test_matrix_stream_SOURCES =		\
	test-matrix-stream.C

test_la_block_lanczos_SOURCES =		\
	test-la-block-lanczos.C

test_mg_block_lanczos_SOURCES =		\
	test-mg-block-lanczos.C

test_rational_solver_adaptive_SOURCES = \
	test-rational-solver-adaptive.C \
	test-common.h

test_smith_form_adaptive_SOURCES=    \
	test-common.h	\
	test-smith-form-adaptive.C

test_fields_SOURCES=    \
	test-common.h	\
	test-fields.C

test_isposdef_SOURCES=    \
	test-common.h	\
	test-isposdef.C

test_ispossemidef_SOURCES=    \
	test-common.h	\
	test-ispossemidef.C
