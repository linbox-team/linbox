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

# All tests should appear in one of the lists of this makefile (and sources given).
# Tests not used for target "check can be in the FULLCHECK_TESTS list.
# All tests will be processed by checker.C (target fullcheck),
# where skips can be specified, warnings issued, package dependency handled.

AM_CPPFLAGS=-I$(top_srcdir) -DLinBoxTestOnly
AM_CXXFLAGS=$(LINBOX_CXXFLAGS) $(FFLAS_FFPACK_CFLAGS) $(LINBOX_DEPS_CFLAGS)
LDADD=$(LINBOX_DEPS_LIBS) $(FFLAS_FFPACK_LIBS) $(LINBOX_LDFLAGS)
AM_LDFLAGS=$(top_builddir)/linbox/liblinbox.la

SUBDIRS = data matrix

PERFPUBLISHERFILE = tests-report.xml

EXTRA_DIST = \
    perfpublisher.sh \
    test-blackbox.h        \
    test-common.h        \
    test-common.inl        \
    test-field.h        \
    test-generic.h        \
    test-matrix-utils.h    \
    test-matrix-domain.h    \
    test-rank.h        \
    test-vector-domain.h    \
    test-smith-form.h       \
    tests.doxy        \
    checker.C

#  tutorial-1: tutorial-1.C

#  test-tutorial$(EXEEXT): test-tutorial.C tutorial-1

# Put an entry for your test-mycomponent.C in one of the test groups.
# Don't forget to add a test_mycomponent_SOURCES entry below.
#
# Once present, tests should remain present, being moved from one group to another.
BASIC_TESTS =            \
    test-solve          \
    test-solve-full             \
    test-smith-form-valence  \
    test-local-smith-form-sparseelim\
    test-smith-form             \
    test-smith-form-adaptive     \
    test-smith-form-iliopoulos  \
    test-smith-form-local        \
    test-last-invariant-factor  \
    test-qlup                    \
    test-det            \
    test-regression        \
    test-regression2       \
    test-rank-ex        \
    test-rank-u32        \
    test-rank-md        \
    test-rank-Int        \
    test-frobenius          \
    test-rational-solver    \
    test-polynomial-matrix\
    test-rational-solver-adaptive \
    test-randiter-nonzero-prime    \
    test-cra            \
    test-blas-matrix        \
    test-charpoly        \
    test-minpoly                \
    test-commentator        \
    test-isposdef        \
    test-ispossemidef       \
    test-givaropoly        \
    test-gf2            \
    test-givaro-zpz        \
    test-givaro-zpzuns        \
    test-givaro-interfaces        \
    test-echelon-form       \
    test-blas-domain-mul        \
    test-blas-domain            \
    test-hadamard-bound     \
    test-fft                    \
    test-serialization

# Really just one or two of these would be enough for target check.
# The rest can be in target fullcheck.
NTL_TESTS =                 \
    test-smith-form         \
    test-ntl-lzz_p          \
    test-ntl-lzz_pe         \
    test-ntl-lzz_px         \
    test-ntl-lzz_pex        \
    test-ntl-zz_p           \
    test-ntl-toeplitz       \
    test-ntl-hankel         \
    test-ntl-sylvester      \
    test-smith-form-iliopoulos  \
    test-frobenius-small    \
    test-toeplitz-det

# All other tests.
# The checker.C determines which of these are built and run in "make fullcheck".
CHECKER_TESTS =                 \
    test-bitonic-sort           \
    test-blackbox-block-container \
    test-block-wiedemann        \
    test-butterfly              \
    test-companion              \
    test-cradomain              \
    test-diagonal               \
    test-dif                    \
    test-direct-sum             \
    test-dyadic-to-rational     \
    test-echelon-form           \
    test-ffpack                 \
    test-fibb                   \
    test-getentry               \
    test-givaropoly             \
    test-gmp-rational           \
    test-hilbert                \
    test-hom                    \
    test-inverse                \
    test-one-invariant-factor   \
    test-matpoly-mult           \
    test-matrix-domain          \
    test-matrix-stream          \
    test-modular                \
    test-modular-balanced-double\
    test-modular-balanced-float \
    test-modular-balanced-int   \
    test-modular-double         \
    test-modular-float          \
    test-modular-int            \
    test-moore-penrose          \
    test-nullspace              \
    test-optimization           \
    test-order-basis            \
    test-param-fuzzy            \
    test-permutation            \
    test-plain-domain           \
    test-polynomial-matrix      \
    test-quad-matrix            \
    test-rational-matrix-factory\
    test-rational-reconstruction-base \
    test-scalar-matrix          \
    test-smith-form-binary      \
    test-solve-nonsingular      \
    test-sparse                 \
    test-subiterator            \
    test-submatrix              \
    test-subvector              \
    test-sum                    \
    test-trace                  \
    test-transpose              \
    test-triplesbb              \
    test-tutorial               \
    test-vector-domain          \
    test-random-matrix          \
    test-zero-one               \
    test-toom-cook              \
    test-toeplitz-det           \
    test-dense

CHECKER_NTL_TESTS =  ${NTL_TESTS}           \
    test-polynomial-ring        \
    test-toeplitz               \
    test-ntl-rr


FULLCHECK_TESTS = ${CHECKER_TESTS} \
    test-weak-popov-form        \
    test-mpi-comm               \
    test-rat-solve              \
    test-rat-minpoly            \
    test-rat-charpoly           \
    test-poly-det               \
    test-sparse-map-map         \
    test-block-ring             \
    test-ftrmm                  \
    test-image-field            \
    test-la-block-lanczos       \
    test-mg-block-lanczos       \
    test-modular-byte           \
    test-modular-short

FULLCHECK_NTL_TESTS = ${CHECKER_NTL_TESTS} \
    test-frobenius-leading-invariants      \
    test-frobenius-large    \
    test-weak-popov-form

# @note Fact is make check won't run the tests with mpirun,
# so it will always fail
# if LINBOX_HAVE_MPI
# MPI_TESTS =     \
#     test-mpi-comm
# endif

if LINBOX_HAVE_NTL
USE_NTL_TESTS =     \
    test-ntl-lzz_pex    \
    test-toeplitz-det
CHECKER_TESTS += ${CHECKER_NTL_TESTS}
FULLCHECK_TESTS += ${FULLCHECK_NTL_TESTS}
endif

# if LINBOX_HAVE_SAGE
# SAGE_LIBS = -L@libdir@ -llinboxsage
# test_regression_LDFLAGS = $(SAGE_LIBS)
# endif

OCL_TESTS = test-opencl-domain

if LINBOX_HAVE_OCL
USE_OCL_TESTS = $(OCL_TESTS)
endif

# check builds and runs these
TESTS =         \
    $(BASIC_TESTS)      \
    $(USE_NTL_TESTS)    \
    $(USE_OCL_TESTS)    \
    $(MPI_TESTS)

# alternate definition of check target
#check: checker
#    checker -c

EXTRA_PROGRAMS =    \
    $(BASIC_TESTS)      \
    $(NTL_TESTS)    \
    $(FULLCHECK_TESTS)  \
    $(OCL_TESTS)

CLEANFILES= checker     \
    $(TESTS)        \
    $(BASIC_TESTS)      \
    $(FULLCHECK_TESTS)      \
    $(NTL_TESTS)        \
    $(OCL_TESTS)        \
    $(PERFPUBLISHERFILE)

test_bitonic_sort_SOURCES =         test-bitonic-sort.C
test_blackbox_block_container_SOURCES = test-blackbox-block-container.C
test_blas_domain_SOURCES =          test-blas-domain.C
test_blas_domain_mul_SOURCES =      test-blas-domain-mul.C
test_blas_matrix_SOURCES =          test-blas-matrix.C
test_block_ring_SOURCES =           test-block-ring.C
test_block_wiedemann_SOURCES =      test-block-wiedemann.C
test_butterfly_SOURCES =        test-butterfly.C test-vector-domain.h test-blackbox.h
test_charpoly_SOURCES =         test-charpoly.C
test_commentator_SOURCES =          test-commentator.C
test_companion_SOURCES =        test-companion.C
test_cradomain_SOURCES =        test-cradomain.C test-common.h
test_cra_SOURCES =              test-cra.C test-common.h
test_dense_SOURCES =            test-dense.C test-common.h
test_det_SOURCES =              test-det.C
test_diagonal_SOURCES =         test-diagonal.C
test_dif_SOURCES =              test-dif.C
test_direct_sum_SOURCES =           test-direct-sum.C
test_dyadic_to_rational_SOURCES =       test-dyadic-to-rational.C
test_echelon_form_SOURCES =         test-echelon-form.C
test_fft_SOURCES =                  test-fft.C
test_ffpack_SOURCES =           test-ffpack.C
test_fibb_SOURCES =             test-fibb.C
test_frobenius_SOURCES =        test-frobenius.C
test_ftrmm_SOURCES =            test-ftrmm.C
test_getentry_SOURCES =         test-getentry.C
test_gf2_SOURCES =              test-gf2.C
test_givaropoly_SOURCES =           test-givaropoly.C
test_givaro_zpz_SOURCES =           test-givaro-zpz.C
test_givaro_zpzuns_SOURCES =        test-givaro-zpzuns.C
test_givaro_interfaces_SOURCES =        test-givaro-interfaces.C
test_gmp_rational_SOURCES =         test-gmp-rational.C
test_hadamard_bound_SOURCES =       test-hadamard-bound.C
test_hilbert_SOURCES =          test-hilbert.C
test_hom_SOURCES =              test-hom.C
test_image_field_SOURCES =          test-image-field.C
test_inverse_SOURCES =          test-inverse.C
test_isposdef_SOURCES =         test-isposdef.C
test_ispossemidef_SOURCES =         test-ispossemidef.C
test_la_block_lanczos_SOURCES =     test-la-block-lanczos.C
test_last_invariant_factor_SOURCES =    test-last-invariant-factor.C
test_matpoly_mult_SOURCES=          test-matpoly-mult.C
test_matrix_domain_SOURCES =        test-matrix-domain.C test-common.h
test_matrix_stream_SOURCES =        test-matrix-stream.C
test_mg_block_lanczos_SOURCES =     test-mg-block-lanczos.C
test_minpoly_SOURCES =          test-minpoly.C
test_modular_balanced_double_SOURCES =  test-modular-balanced-double.C
test_modular_balanced_float_SOURCES =   test-modular-balanced-float.C
test_modular_balanced_int_SOURCES =     test-modular-balanced-int.C
test_modular_byte_SOURCES =             test-modular-byte.C
test_modular_double_SOURCES =           test-modular-double.C
test_modular_float_SOURCES =            test-modular-float.C
test_modular_int_SOURCES =              test-modular-int.C
test_modular_short_SOURCES =            test-modular-short.C
test_modular_SOURCES =                  test-modular.C
test_moore_penrose_SOURCES =            test-moore-penrose.C
test_ntl_hankel_SOURCES =               test-ntl-hankel.C
test_ntl_lzz_pe_SOURCES =               test-ntl-lzz_pe.C test-field.h
test_ntl_lzz_pex_SOURCES =              test-ntl-lzz_pex.C test-field.h
test_ntl_lzz_p_SOURCES =                test-ntl-lzz_p.C test-field.h
test_ntl_lzz_px_SOURCES =               test-ntl-lzz_px.C test-field.h
test_ntl_rr_SOURCES =                   test-ntl-rr.C
test_ntl_sylvester_SOURCES =            test-ntl-sylvester.C
test_ntl_toeplitz_SOURCES =             test-ntl-toeplitz.C
test_ntl_zz_p_SOURCES =                 test-ntl-zz_p.C
test_nullspace_SOURCES =                test-nullspace.C
test_one_invariant_factor_SOURCES =     test-one-invariant-factor.C
test_opencl_domain_SOURCES =            test-opencl-domain.C
test_optimization_SOURCES =             test-optimization.C
test_order_basis_SOURCES =              test-order-basis.C
test_param_fuzzy_SOURCES =              test-param-fuzzy.C
test_permutation_SOURCES =              test-permutation.C
test_polynomial_matrix_SOURCES=         test-polynomial-matrix.C
test_plain_domain_SOURCES =             test-plain-domain.C
test_poly_det_SOURCES =                 test-poly-det.C
test_qlup_SOURCES =                     test-qlup.C
test_quad_matrix_SOURCES =              test-quad-matrix.C
test_randiter_nonzero_prime_SOURCES =   test-randiter-nonzero-prime.C
test_random_matrix_SOURCES =        test-random-matrix.C
test_rank_ex_SOURCES =         test-rank-ex.C
test_rank_Int_SOURCES =         test-rank-Int.C test-rank.h
test_rank_md_SOURCES =          test-rank-md.C test-rank.h
test_rank_u32_SOURCES =         test-rank-u32.C test-rank.h
test_rat_charpoly_SOURCES =         test-rat-charpoly.C test-common.h
test_rational_matrix_factory_SOURCES =  test-rational-matrix-factory.C
test_rational_reconstruction_base_SOURCES = test-rational-reconstruction-base.C
test_rational_solver_adaptive_SOURCES = test-rational-solver-adaptive.C test-common.h
test_rational_solver_SOURCES =      test-rational-solver.C
test_rat_minpoly_SOURCES =          test-rat-minpoly.C test-common.h
test_rat_solve_SOURCES =        test-rat-solve.C test-common.h
test_regression_SOURCES =           test-regression.C
test_regression2_SOURCES =           test-regression2.C
test_scalar_matrix_SOURCES =        test-scalar-matrix.C
test_serialization_SOURCES =         test-serialization.C
test_smith_form_adaptive_SOURCES =      test-smith-form-adaptive.C test-common.h
test_smith_form_binary_SOURCES =    test-smith-form-binary.C
test_smith_form_iliopoulos_SOURCES =    test-smith-form-iliopoulos.C
test_smith_form_local_SOURCES =     test-smith-form-local.C
test_smith_form_valence_SOURCES = test-smith-form-valence.C
test_local_smith_form_sparseelim_SOURCES = test-local-smith-form-sparseelim.C
test_smith_form_SOURCES =           test-smith-form.C
test_solve_nonsingular_SOURCES =    test-solve-nonsingular.C
test_solve_SOURCES =            test-solve.C
test_solve_full_SOURCES =               test-solve-full.C
test_sparse_SOURCES =           test-sparse.C test-common.h
test_sparse_map_map_SOURCES =         test-sparse-map-map.C test-blackbox.h
test_subiterator_SOURCES =          test-subiterator.C test-common.h
test_submatrix_SOURCES =        test-submatrix.C test-common.h
test_subvector_SOURCES =        test-subvector.C test-common.h
test_sum_SOURCES =              test-sum.C
test_toeplitz_det_SOURCES =         test-toeplitz-det.C
test_toom_cook_SOURCES =        test-toom-cook.C
test_trace_SOURCES =            test-trace.C
test_transpose_SOURCES =        test-transpose.C
#test_triplesbb_omp_SOURCES =        test-triplesbb-omp.C
test_triplesbb_SOURCES =        test-triplesbb.C
test_tutorial_SOURCES =         test-tutorial.C
test_vector_domain_SOURCES =        test-vector-domain.C test-vector-domain.h
test_zero_one_SOURCES =         test-zero-one.C
test_polynomial_ring_SOURCES =      test-polynomial-ring.C
test_frobenius_leading_invariants_SOURCES =    test-frobenius-leading-invariants.C
test_frobenius_small_SOURCES =      test-frobenius-small.C
test_frobenius_large_SOURCES =      test-frobenius-large.C
test_weak_popov_form_SOURCES =      test-weak-popov-form.C
test_mpi_comm_SOURCES =         test-mpi-comm.C
test_toeplitz_SOURCES =                 test-toeplitz.C
checker_SOURCES      =    checker.C 

# Perfpublisher script interaction - AB 2014/12/11
perfpublisher:
	$(top_srcdir)/tests/perfpublisher.sh "$(PERFPUBLISHERFILE)" "$(TESTS)" "$(CXX)"
	mv $(PERFPUBLISHERFILE) ..
	make clean
	mv ../$(PERFPUBLISHERFILE) .

# for compilation of non-enumerated programs
%:%.C
	$(CXX) $(CXXFLAGS) $(AM_CPPFLAGS)  ${INCLUDES} $(AM_CXXFLAGS) $*.C -o $@ $(LDADD)

fullcheck: ${TESTS} ${CHECKER_TESTS} checker
	./checker
