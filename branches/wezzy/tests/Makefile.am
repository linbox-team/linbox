# Copyright (c) 2010-2011 the LinBox group
# ========LICENCE========
# This file is part of the library LinBox.
#
# LinBox is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========




AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/linbox
# AM_LDFLAGS = -llapack -lcblas -latlas
#  AM_LDFLAGS = $(top_srcdir)/linbox/.libs/liblinbox.a
AM_LDFLAGS = $(top_srcdir)/linbox/liblinbox.la
#  SAFER_FLAGS = -g -Wall -Wextra -Wno-unused-parameter
# bds: no need to optimize for unit/regression tests, and opt adds a lot to compilation time
#AM_CXXFLAGS = @TESTS_CFLAGS@ -DLinBoxSrcOnly
#  AM_CXXFLAGS = $(SAFER_FLAGS) -DLinBoxSrcOnly
#AM_CXXFLAGS = $(SAFER_FLAGS) -DLinBoxSrcOnly -O2
#  AM_CXXFLAGS = -g -DLinBoxSrcOnly -Wall -DDEBUG -O0
AM_CXXFLAGS = -DLinBoxSrcOnly  @TESTS_CFLAGS@ @DEBUG_CFLAGS@
# AM_CXXFLAGS = -g -DLinBoxSrcOnly -Wall -D__LINBOX_HAVE_DGETRF -D__LINBOX_HAVE_DGETRI
# AM_CFLAGS = -g -Wall

SUBDIRS = data matrix

AM_CPPFLAGS+= $(CBLAS_FLAG) $(FFLAS_FFPACK_CFLAGS) $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(IML_CFLAGS) $(FLINT_CFLAGS) $(OCL_CFLAGS)
LDADD = $(GIVARO_LIBS) $(BLAS_LIBS) $(NTL_LIBS)  $(GMP_LIBS) $(LIDIA_LIBS) $(MPFR_LIBS) $(FPLLL_LIBS)  $(IML_LIBS) $(FLINT_LIBS) $(LDFLAGS) $(OCL_LIBS)
#  AM_LDFLAGS+=-static -Wl,-rpath -Wl,$(GIVARO_LIBS)

EXTRA_DIST=            \
		test-common.C  \
		tests.doxy     \
		checker.C

#  test-generic-for-quad.h


tutorial-1: tutorial-1.C

test-tutorial$(EXEEXT): test-tutorial.C tutorial-1

# Put an entry for your test-mycomponent.C here (or in one of the other test groups).
# Don't forget to add a test_mycomponent_SOURCES entry below.
#
# Once present, tests should remain present, being moved from one group to another.
# Note that among the available groups are NONCOMPILING_TESTS and FAILING_TESTS, which
# will be ignored by "make check" and can be temporary storage for problematic tests.
BASIC_TESTS =  \
        test-commentator              \
        test-det                      \
        test-frobenius                \
        test-rank                     \
        test-solve                    \
        test-nullspace                \
        test-smith-form-binary        \
        test-rational-solver          \
        test-rational-solver-adaptive \
        test-charpoly                 \
        test-cra

FULLCHECK_TESTS =                 \
        test-bitonic-sort             \
        test-blackbox-block-container \
        test-blas-domain              \
		test-block-wiedemann		  \
        test-butterfly                \
        test-companion                \
        test-dense                    \
        test-diagonal                 \
        test-dif                      \
        test-direct-sum               \
        test-echelon-form             \
        test-ffpack                   \
        test-getentry                 \
        test-gf2                      \
        test-gmp-rational             \
        test-hilbert                  \
        test-hom                      \
        test-inverse                  \
        test-isposdef                 \
        test-ispossemidef             \
        test-matrix-domain            \
        test-matrix-stream            \
        test-minpoly                  \
        test-modular                  \
        test-modular-balanced-int     \
        test-modular-balanced-float   \
        test-modular-balanced-double  \
        test-modular-byte             \
        test-modular-double           \
        test-modular-float            \
        test-modular-int              \
        test-modular-short            \
        test-block-ring               \
        test-PID-integer              \
        test-moore-penrose            \
        test-randiter-nonzero         \
        test-rational-matrix-factory  \
        test-scalar-matrix            \
        test-solve-nonsingular        \
        test-sparse                   \
        test-subiterator              \
        test-submatrix                \
        test-subvector                \
        test-sum                      \
        test-trace                    \
        test-transpose              \
        test-triplesbb                \
        test-tutorial                \
        test-vector-domain            \
        test-qlup                     \
        test-zero-one                 \
        test-last-invariant-factor    \
        test-cradomain                \
        test-rational-reconstruction-base

#  test-integer-tools

# tests of components not adhering to LinBox requirements.
NONCOMPLIANT_TESTS =                   \
        test-ntl-RR                    \
        test-param-fuzzy

#These need clapack
#test-rational-solver-adaptive

DUBIOUS_TESTS =                 \
	test-unparametric-field       \
	test-bmseq		\
		test-zo
#  test-quad-matrix

NONCOMPILING_TESTS =           \
        test-dense-zero-one

FAILING_TESTS =                  \
        test-rat-charpoly        \
		test-rat-solve           \
		test-rat-minpoly         \
		test-smith-form-local    \
		test-smith-form-adaptive \
		test-dyadic-to-rational  \
		test-mg-block-lanczos    \
		test-ftrmm               \
        test-la-block-lanczos

USE_NONCOMPILING_TESTS =
USE_FAILING_TESTS =
# comment out these two lines to avoid facing reality
#USE_NONCOMPILING_TESTS += $(NONCOMPILING_TESTS)
#USE_FAILING_TESTS += $(FAILING_TESTS)

# a benchmarker, not to be included in check.
BENCHMARKS =                    \
        benchmark-fields        \
        benchmark-blas-domain

# historical tests, no longer in use
#        test-nag-sparse                  replaced by tripplesbb
#       test-lsp         lsp.h superceded by factorized-matrix.h
# removed test
#        test-nag-sparse
#        test-image-field       replaced by hom
#       test-lsp
# I put test-gssv_rank in hmrg:~saunders/gen_superlu.  It is not linbox until and unless it is made to work
# test-gssv is an ntl test


# no explicit test needed, I guess...
#        Transpose is tested in test-triplesbb
#        Compose is tested in test-butterfly

GIVARO_TESTS =                        \
        test-givaro-fields

NTL_TESTS =                           \
        test-smith-form               \
        test-ntl-lzz_p                \
        test-ntl-ZZ_p                 \
        test-ntl-toeplitz             \
        test-ntl-hankel               \
        test-ntl-sylvester            \
        test-smith-form-iliopoulos    \
        test-toeplitz-det

LIDIA_TESTS =                         \
        test-lidia-gfq

ATLAS_TESTS =                         \
        test-optimization

BLAS3_UNIT_TESTS = \
		   test-toom-cook

OCL_TESTS =                           \
        test-opencl-domain

if LINBOX_HAVE_NTL
USE_NTL_TESTS = $(NTL_TESTS)
endif

if LINBOX_HAVE_GIVARO
USE_GIVARO_TESTS = $(GIVARO_TESTS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_TESTS = $(LIDIA_TESTS)
endif

if LINBOX_HAVE_OCL
USE_OCL_TESTS = $(OCL_TESTS)
endif

USE_ATLAS_TESTS = $(ATLAS_TESTS)


TESTS =                               \
        $(BASIC_TESTS)                \
        $(USE_GIVARO_TESTS)           \
        $(USE_NTL_TESTS)              \
        $(USE_LIDIA_TESTS)            \
        $(USE_NONCOMPILING_TESTS)     \
        $(USE_FAILING_TESTS)          \
        $(USE_ATLAS_TESTS)            \
        $(USE_OCL_TESTS)



EXTRA_PROGRAMS =                      \
        $(BASIC_TESTS)                \
        $(GIVARO_TESTS)               \
        $(NTL_TESTS)                  \
        $(ATLAS_TESTS)                \
        $(LIDIA_TESTS)                \
        $(BENCHMARKS)                 \
        $(NONEXISTING_TESTS)          \
        $(NONCOMPILING_TESTS)         \
        $(NONCOMPLIANT_TESTS)         \
        $(DUBIOUS_TESTS)              \
        $(FULLCHECK_TESTS)            \
        $(FAILING_TESTS)              \
        $(BLAS3_UNIT_TESTS)           \
        $(OCL_TESTS)


test_commentator_SOURCES =                \
        test-commentator.C                \
        test-common.h

test_modular_SOURCES =                        \
        test-modular.C                        \
        test-generic.h                        \
        test-common.h

test_blackbox_block_container_SOURCES =   \
	test-blackbox-block-container.C    \
	test-common.h

test_solve_nonsingular_SOURCES =   \
	test-solve-nonsingular.C   \
	test-common.h

test_block_ring_SOURCES =                 \
        test-block-ring.C                 \
        test-field.h                      \
        test-common.h

test_PID_integer_SOURCES =                 \
        test-PID-integer.C                 \
        test-field.h                      \
        test-common.h

test_hom_SOURCES =                        \
        test-hom.C                        \
        test-generic.h                    \
        test-common.h

test_gf2_SOURCES =                        \
        test-gf2.C                        \
        test-generic.h                    \
        test-common.h

test_gmp_rational_SOURCES =               \
        test-gmp-rational.C               \
        test-generic.h                    \
        test-common.h

test_unparametric_field_SOURCES =        \
        test-unparametric-field.C        \
        test-generic.h                   \
        test-common.h

test_randiter_nonzero_SOURCES =          \
        test-randiter-nonzero.C          \
        test-common.h

test_rational_matrix_factory_SOURCES =  \
        test-rational-matrix-factory.C  \
        test-common.h

test_rational_reconstruction_base_SOURCES =  \
	  test-rational-reconstruction-base.C  \
	  test-common.h

test_rat_solve_SOURCES =                \
        test-rat-solve.C                \
        test-common.h

test_rat_charpoly_SOURCES =             \
        test-rat-charpoly.C             \
        test-common.h

test_rat_minpoly_SOURCES =              \
        test-rat-minpoly.C              \
        test-common.h

test_subiterator_SOURCES =              \
        test-subiterator.C              \
        test-common.h

test_subvector_SOURCES =                \
        test-subvector.C                \
        test-common.h

test_vector_domain_SOURCES =            \
        test-vector-domain.C            \
        test-vector-domain.h

test_matrix_domain_SOURCES =            \
        test-matrix-domain.C            \
        test-common.h

test_dense_SOURCES =                    \
        test-dense.C                    \
        test-common.h

test_sparse_SOURCES =                   \
        test-sparse.C                   \
        test-common.h

test_submatrix_SOURCES =                \
        test-submatrix.C                \
        test-common.h

test_scalar_matrix_SOURCES =            \
        test-scalar-matrix.C            \
        test-generic.h                  \
        test-common.h

test_bmseq_SOURCES =            \
        test-bmseq.C            \
        test-common.h

test_butterfly_SOURCES =                \
        test-vector-domain.h            \
        test-blackbox.h                 \
        test-butterfly.C                \
        test-generic.h                  \
        test-common.h

test_charpoly_SOURCES =                 \
        test-charpoly.C                 \
        test-common.h

test_diagonal_SOURCES =                 \
        test-diagonal.C                 \
        test-generic.h                  \
        test-common.h

test_sum_SOURCES =                      \
        test-sum.C                      \
        test-generic.h                  \
        test-common.h

test_dif_SOURCES =                          \
        test-dif.C        test-generic.h    \
        test-common.h

test_triplesbb_SOURCES =                \
        test-triplesbb.C                \
        test-generic.h                  \
        test-common.h

test_tutorial_SOURCES =                \
        test-tutorial.C

test_hilbert_SOURCES =                        \
        test-hilbert.C        test-generic.h        \
        test-common.h

test_ntl_toeplitz_SOURCES =                \
        test-ntl-toeplitz.C                \
        test-common.h

test_ntl_hankel_SOURCES =                \
        test-ntl-hankel.C                \
        test-common.h


test_ntl_sylvester_SOURCES =                \
        test-ntl-sylvester.C                \
        test-common.h

test_toeplitz_det_SOURCES =                \
        test-toeplitz-det.C

test_trace_SOURCES =                         \
        test-trace.C                         \
        test-common.h

test_getentry_SOURCES =                         \
        test-getentry.C                         \
        test-common.h

test_minpoly_SOURCES =                         \
        test-minpoly.C                         \
        test-common.h

test_det_SOURCES =                        \
        test-det.C                        \
        test-common.h

test_rank_SOURCES =                        \
        test-rank.C                        \
        test-common.h

test_qlup_SOURCES =                        \
        test-qlup.C                        \
        test-common.h

test_solve_SOURCES =                        \
        test-solve.C                        \
        test-common.h

test_direct_sum_SOURCES =                \
        test-direct-sum.C                \
        test-generic.h                   \
        test-common.h

test_echelon_form_SOURCES =                \
        test-echelon-form.C                \
        test-common.h

test_frobenius_SOURCES =                \
        test-frobenius.C                \
        test-generic.h                  \
        test-common.h

test_companion_SOURCES =                \
        test-companion.C                \
        test-generic.h                  \
        test-common.h

test_inverse_SOURCES =                        \
        test-inverse.C                        \
        test-common.h

test_moore_penrose_SOURCES =                \
        test-moore-penrose.C                \
        test-common.h

test_ntl_RR_SOURCES =                        \
        test-ntl-RR.C                        \
        test-generic.h                       \
        test-common.h

test_ntl_ZZ_p_SOURCES =                        \
        test-ntl-ZZ_p.C                        \
        test-generic.h                         \
        test-common.h

test_ntl_lzz_p_SOURCES =                \
        test-ntl-lzz_p.C                \
        test-generic.h                  \
        test-common.h

test_givaro_fields_SOURCES =                \
        test-givaro-fields.C                \
        test-generic.h                   \
        test-common.h

test_lidia_gfq_SOURCES =                \
        test-lidia-gfq.C                \
        test-generic.h                  \
        test-common.h

test_param_fuzzy_SOURCES =                \
        test-param-fuzzy.C                \
        test-generic.h                    \
        test-common.h

#  test_quad_matrix_SOURCES =                \
#  test-quad-matrix.C                \
#  test-generic-for-quad.h           \
#  test-generic.h                    \
#  test-common.h

test_zo_SOURCES =                        \
        test-zo.C                        \
        test-generic.h                   \
        test-common.h

test_zero_one_SOURCES =                       \
        test-zero-one.C                       \
        test-generic.h                        \
        test-common.h

test_modular_int_SOURCES =                \
        test-modular-int.C                \
        test-common.h

test_modular_short_SOURCES =                \
        test-modular-short.C                \
        test-generic.h                      \
        test-common.h

test_modular_double_SOURCES =                \
        test-modular-double.C                \
        test-generic.h                       \
        test-common.h

test_modular_float_SOURCES =                \
        test-modular-float.C                \
        test-generic.h                       \
        test-common.h


test_modular_byte_SOURCES =          \
        test-modular-byte.C          \
        test-generic.h               \
        test-common.h

test_modular_balanced_int_SOURCES =  \
        test-modular-balanced-int.C  \
        test-generic.h               \
        test-common.h

test_modular_balanced_double_SOURCES =  \
        test-modular-balanced-double.C  \
        test-generic.h               \
        test-common.h

test_modular_balanced_float_SOURCES =  \
        test-modular-balanced-float.C  \
        test-generic.h               \
        test-common.h

test_smith_form_iliopoulos_SOURCES=  \
        test-common.h                \
        test-smith-form-iliopoulos.C

test_smith_form_binary_SOURCES=  \
        test-common.h            \
        test-smith-form-binary.C

test_smith_form_SOURCES=    \
        test-common.h       \
        test-smith-form.C

test_last_invariant_factor_SOURCES=    \
        test-common.h                  \
        test-last-invariant-factor.C

test_smith_form_local_SOURCES=    \
        test-common.h             \
        test-smith-form-local.C

test_optimization_SOURCES=    \
        test-common.h         \
        test-optimization.C

test_nullspace_SOURCES=    \
        test-common.h      \
        test-nullspace.C

test_bitonic_sort_SOURCES=              \
        test-common.h                   \
        test-bitonic-sort.C

test_rational_solver_SOURCES =          \
        test-common.h                   \
        test-rational-solver.C

test_blas_domain_SOURCES=               \
        test-common.h                   \
        test-blas-domain.C

test_ffpack_SOURCES=                    \
        test-common.h                   \
        test-ffpack.C

test_matrix_stream_SOURCES =            \
        test-matrix-stream.C

test_la_block_lanczos_SOURCES =         \
        test-la-block-lanczos.C

test_mg_block_lanczos_SOURCES =         \
        test-mg-block-lanczos.C

test_block_wiedemann_SOURCES =         \
        test-block-wiedemann.C

test_rational_solver_adaptive_SOURCES = \
        test-rational-solver-adaptive.C \
        test-common.h

test_dyadic_to_rational_SOURCES = \
        test-dyadic-to-rational.C \
        test-common.h

test_cra_SOURCES =                      \
        test-cra.C                      \
        test-common.h

test_cradomain_SOURCES =                \
        test-cradomain.C                \
        test-common.h

# NTL test

test_smith_form_adaptive_SOURCES=      \
        test-common.h                  \
        test-smith-form-adaptive.C

test_transpose_SOURCES =                \
        test-transpose.C                \
        test-common.h

benchmark_fields_SOURCES=              \
        test-common.h                  \
        benchmark-fields.C

benchmark_blas_domain_SOURCES=         \
        test-common.h                  \
        benchmark-blas-domain.C

test_isposdef_SOURCES=                 \
        test-common.h                  \
        test-isposdef.C

test_dense_zero_one_SOURCES=           \
        test-blackbox.h                \
        test-dense-zero-one.C

test_ispossemidef_SOURCES=            \
        test-common.h                 \
        test-ispossemidef.C

test_ftrmm_SOURCES =                  \
        test-ftrmm.C

# OCL Test

test_opencl_domain_SOURCES=           \
        test-common.h                 \
        test-opencl-domain.C

test_toom_cook_SOURCES = test-toom-cook.C

#  test_integer_tools_SOURCES =          \
#  test-integer-tools.C

#check: checker
#	checker -c

# for compilation of new non-enumerated tests
%:%.C
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDADD) -lOpenCL

fullcheck: checker
	./checker

checker: checker.C

CLEANFILES= checker               \
			$(TESTS)              \
			$(BASIC_TESTS)        \
			$(FAILING_TESTS)      \
			$(NONCOMPLIANT_TESTS) \
			$(BENCHMARKS)         \
			$(DUBIOUS_TESTS)      \
			$(FULLCHECK_TESTS)      \
			$(NTL_TESTS)          \
			$(GIVARO_TESTS)       \
			$(ATLAS_TESTS)

