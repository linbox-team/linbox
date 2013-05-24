# Copyright (c) 2010 the LinBox group
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



pkgincludesubdir=$(pkgincludedir)/algorithms

SUBDIRS=gauss matrix-blas3 opencl-kernels IML

noinst_LTLIBRARIES=libalgorithms.la

libalgorithms_la_SOURCES= diophantine-solver.C

#  AM_CPPFLAGS= $(CBLAS_FLAG) $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)

AM_CXXFLAGS = @DEFAULT_CFLAGS@ -DDISABLE_COMMENTATOR $(CBLAS_FLAG) $(FFLAS_FFPACK_CFLAGS) $(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(FPLLL_CFLAGS) $(OCL_CFLAGS)

OCL_HDRS =                                  \
    opencl-domain.h                         \
    opencl-domain.inl                       \
    opencl-domain-util.inl                  \
    opencl-domain-memory.inl                \
    opencl-domain-factory.h                 \
    opencl-resource-controller.h            \
    opencl-environ.h

if LINBOX_HAVE_OCL
USE_OCL_HDRS = $(OCL_HEADERS)
endif

pkgincludesub_HEADERS =                     \
    blackbox-container.h                    \
    blackbox-container-symmetric.h          \
    blackbox-container-symmetrize.h         \
    blackbox-container-base.h               \
    bbsolve.h                               \
    double-det.h                            \
    det-rational.h                          \
    massey-domain.h                         \
    wiedemann.h wiedemann.inl               \
    lanczos.h   lanczos.inl                 \
    block-lanczos.h block-lanczos.inl       \
    mg-block-lanczos.h mg-block-lanczos.inl \
    la-block-lanczos.h la-block-lanczos.inl \
    eliminator.h    eliminator.inl          \
    gauss.h                                 \
    iml.h                                   \
    triangular-solve.h                      \
    gauss-gf2.h                             \
    triangular-solve-gf2.h                  \
    dense-container.h                       \
    cra-mpi.h                               \
    cra-kaapi.h                             \
    cra-domain.h                            \
    cra-domain-seq.h                        \
    cra-domain-omp.h                        \
    cra-early-multip.h                      \
    cra-early-single.h                      \
    cra-full-multip.h                       \
    cra-full-multip-fixed.h                 \
    cra-givrnsfixed.h                       \
    lazy-product.h                          \
    rational-cra.h                          \
    rational-cra2.h                         \
    rational-cra-early-multip.h             \
    rational-cra-early-single.h             \
    rational-cra-full-multip.h              \
    blas-domain.h                           \
    blas-domain.inl                         \
    blas-domain-multimod.h                  \
    bitonic-sort.h                          \
    blackbox-block-container-base.h         \
    blackbox-block-container.h              \
    block-massey-domain.h                   \
    block-wiedemann.h                       \
    default.h                               \
    signature.h                             \
    smith-form-iliopoulos.h                 \
    one-invariant-factor.h                  \
    last-invariant-factor.h                 \
    hybrid-det.h                            \
    lifting-container.h                     \
    smith-form-local.h                      \
    smith-form-local2.inl                   \
    matrix-inverse.h                        \
    matrix-hom.h                            \
    matrix-rank.h                           \
    numeric-solver-lapack.h                 \
    rational-solver-sn.h                    \
    rational-solver-sn.inl                  \
    dyadic-to-rational.h                    \
    rational-reconstruction.h               \
    rational-reconstruction-base.h          \
    classic-rational-reconstruction.h       \
    fast-rational-reconstruction.h          \
    rational-solver.h                       \
    rational-solver2.h                      \
    rational-solver.inl                     \
    diophantine-solver.h                    \
    diophantine-solver.inl                  \
    smith-form-binary.h                     \
    smith-form-adaptive.h                   \
    smith-form-adaptive.inl                 \
    smith-form-sparseelim-local.h           \
    smith-form-sparseelim-poweroftwo.h      \
    rational-reconstruction2.h  \
    rational-solver-adaptive.h  \
    varprec-cra-early-single.h  \
    varprec-cra-early-multip.h  \
    vector-fraction.h           \
    whisart_trace.h             \
    cia.h                       \
    lattice.h lattice.inl      \
    bbcharpoly.h                \
    charpoly-rational.h         \
    minpoly-integer.h           \
    minpoly-rational.h          \
    sigma-basis.h               \
    matpoly-mult.h              \
    echelon-form.h              \
    toeplitz-det.h              \
    dense-nullspace.h dense-nullspace.inl   \
    linbox-tags.h                           \
    short-vector.h                          \
    rns.h rns.inl                           \
    $(USE_OCL_HDRS)

EXTRA_DIST =           \
    algorithms.doxy    \
    cra.doxy           \
    iml.doxy           \
    wiedemann.doxy     \
    lanczos.doxy       \
    sigmabasis.doxy    \
    elimination.doxy   \
    lattice.doxy       \
    diophantine.doxy   \
    gauss.doxy


