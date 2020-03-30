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

AM_CPPFLAGS=-I$(top_srcdir)
AM_CXXFLAGS=$(LINBOX_CXXFLAGS) $(FFLAS_FFPACK_CFLAGS) $(LINBOX_DEPS_CFLAGS)
LDADD=$(LINBOX_DEPS_LIBS) $(FFLAS_FFPACK_LIBS) $(LINBOX_LDFLAGS)

SUBDIRS=dixon-solver gauss opencl-kernels matrix-blas3 polynomial-matrix
#  IML

noinst_LTLIBRARIES=libalgorithms.la

libalgorithms_la_SOURCES= diophantine-solver.C

OCL_HDRS =                                  \
    opencl-resource-controller.h            \
    opencl-environ.h

#  if LINBOX_HAVE_OCL // we include them anyway. (for make dist when OCL not available)
USE_OCL_HDRS = $(OCL_HDRS)
#  endif



pkgincludesub_HEADERS =            \
	alt-blackbox-block-container.h     \
	bbcharpoly.h                       \
	bitonic-sort.h                     \
	blackbox-block-container-base.h    \
	blackbox-block-container.h         \
	blackbox-container-base.h          \
	blackbox-container.h               \
	blackbox-container-symmetric.h     \
	blackbox-container-symmetrize.h    \
	block-coppersmith-domain.h         \
	block-lanczos.h                    \
	block-lanczos.inl                  \
	block-massey-domain.h              \
	block-wiedemann.h                  \
	charpoly-rational.h                \
	cia.h                              \
	classic-rational-reconstruction.h  \
	coppersmith.h                      \
	coppersmith-invariant-factors.h    \
	cra-domain.h                       \
	cra-domain-omp.h                   \
	cra-domain-sequential.h                   \
	cra-builder-early-multip.h                 \
	cra-builder-full-multip-fixed.h            \
	cra-builder-full-multip.h                  \
	cra-givrnsfixed.h                  \
	cra-kaapi.h                        \
	cra-distributed.h                  \
	cra-builder-single.h                       \
	default.h                          \
	dense-container.h                  \
	dense-nullspace.h                  \
	dense-nullspace.inl                \
	det-rational.h                     \
	diophantine-solver.h               \
	diophantine-solver.inl             \
	double-det.h                       \
	dyadic-to-rational.h               \
	eliminator.h                       \
	eliminator.inl                     \
	fast-rational-reconstruction.h     \
	frobenius-large.h                  \
	frobenius-small.h                  \
	gauss-gf2.h                        \
	gauss.h                            \
	hybrid-det.h                       \
	invariant-factors.h                \
	invert-tb.h                        \
	la-block-lanczos.h                 \
	la-block-lanczos.inl               \
	lanczos.h                          \
	lanczos.inl                        \
	last-invariant-factor.h            \
	lattice.h                          \
	lattice.inl                        \
	lazy-product.h                     \
	lifting-container.h                \
	massey-domain.h                    \
	matpoly-mult.h                     \
	matrix-hom.h                       \
	matrix-inverse.h                   \
	mg-block-lanczos.h                 \
	mg-block-lanczos.inl               \
	minpoly-integer.h                  \
	minpoly-rational.h                 \
	numeric-solver-lapack.h            \
	one-invariant-factor.h             \
	poly-det.h                         \
	poly-dixon.h                       \
	poly-interpolation.h               \
	poly-smith-form.h                  \
	poly-smith-form-local-x.h          \
	rational-cra-var-prec.h                    \
	rational-cra-builder-early-multip.h        \
	rational-cra-builder-early-single.h        \
	rational-cra-builder-full-multip.h         \
	rational-cra.h                     \
	rational-reconstruction2.h         \
	rational-reconstruction-base.h     \
	rational-reconstruction.h          \
	rational-solver-adaptive.h         \
	rational-solver.h                  \
	rational-solver.inl                \
	rational-solver-sn.h               \
	rational-solver-sn.inl             \
	rns.h                              \
	rns.inl                            \
	short-vector.h                     \
	sigma-basis.h                      \
	signature.h                        \
	smith-form-adaptive.h              \
	smith-form-adaptive.inl            \
	smith-form-binary.h                \
	smith-form-iliopoulos.h            \
	smith-form-kannan-bachem.h         \
	smith-form-local2.inl              \
	smith-form-local.h                 \
	smith-form-valence.h               \
	smith-form-sparseelim-local.h      \
	smith-form-sparseelim-poweroftwo.h \
	toeplitz-det.h                     \
	triangular-solve-gf2.h             \
	triangular-solve.h                 \
	cra-builder-var-prec-early-multip.h         \
	cra-builder-var-prec-early-single.h         \
	vector-fraction.h                  \
	weak-popov-form.h                  \
	whisart_trace.h                    \
	wiedemann.h                        \
	wiedemann.inl                      \
	$(USE_OCL_HDRS)

#  iml.h                              \
#

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


