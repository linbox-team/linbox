# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/algorithms

noinst_LTLIBRARIES=libalgorithms.la

libalgorithms_la_SOURCES= diophantine-solver.C

AM_CPPFLAGS= $(CBLAS_FLAG) $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)

pkgincludesub_HEADERS = \
	blackbox-container.h			\
	blackbox-container-symmetric.h		\
	blackbox-container-symmetrize.h		\
	blackbox-container-base.h		\
	bbsolve.h				\
	double-det.h				\
	det-rational.h				\
	massey-domain.h				\
	wiedemann.h	wiedemann.inl		\
	lanczos.h	lanczos.inl		\
	block-lanczos.h	block-lanczos.inl	\
	mg-block-lanczos.h mg-block-lanczos.inl	\
	la-block-lanczos.h la-block-lanczos.inl	\
	eliminator.h	eliminator.inl		\
	gauss.h					\
	gauss.inl				\
	gauss-det.inl				\
	gauss-rank.inl				\
	gauss-solve.inl				\
	gauss-nullspace.inl			\
	gauss-elim.inl				\
	gauss-pivot.inl				\
	triangular-solve.h			\
	gauss-gf2.h				\
	gauss-gf2.inl				\
	gauss-elim-gf2.inl			\
	gauss-rank-gf2.inl			\
        gauss-pivot-gf2.inl			\
        gauss-solve-gf2.inl			\
	triangular-solve-gf2.h			\
	dense-container.h			\
	cra-mpi.h				\
	cra-kaapi.h				\
	cra-domain.h				\
	cra-domain-seq.h				\
	cra-domain-omp.h				\
	cra-early-multip.h			\
	cra-early-single.h			\
	cra-full-multip.h			\
	lazy-product.h				\
	rational-cra.h				\
	rational-cra2.h				\
	rational-cra-early-multip.h		\
	rational-cra-early-single.h		\
	rational-cra-full-multip.h		\
	blas-domain.h				\
	blas-domain.inl				\
	bitonic-sort.h				\
	blackbox-block-container-base.h		\
	blackbox-block-container.h		\
	block-massey-domain.h			\
	block-wiedemann.h			\
	default.h				\
	signature.h				\
	smith-form-iliopoulos.h		\
	one-invariant-factor.h			\
	last-invariant-factor.h			\
	hybrid-det.h				\
	lifting-container.h			\
	smith-form-local.h				\
	smith-form-local2.h				\
	matrix-inverse.h			\
	matrix-hom.h				\
	matrix-rank.h				\
	rational-reconstruction.h		\
	rational-reconstruction-base.h		\
	classic-rational-reconstruction.h	\
	rational-solver.h			\
	rational-solver2.h \
	rational-solver.inl			\
	diophantine-solver.h			\
	diophantine-solver.inl			\
	smith-form-binary.h				\
	smith-form-adaptive.h			\
	smith-form-adaptive.inl			\
	smith-form-sparseelim-local.h		\
	rational-reconstruction2.h		\
	rational-solver-adaptive.h		\
	varprec-cra-early-single.h		\
	varprec-cra-early-multip.h		\
	vector-fraction.h	\
	whisart_trace.h		\
	cia.h			\
	bbcharpoly.h		\
	charpoly-rational.h	\
	minpoly-integer.h	\
	minpoly-rational.h	\
	sigma-basis.h		\
	matpoly-mult.h		\
	echelon-form.h		\
	toeplitz-det.h 

EXTRA_DIST = algorithms.doxy

