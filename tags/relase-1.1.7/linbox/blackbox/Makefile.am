# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/blackbox

BASIC_HDRS =			\
	archetype.h		\
	factory.h		\
	butterfly.h		\
	hilbert.h		\
	compose.h		\
	permutation.h		\
	squarize.h		\
	sparse1.h		\
	sparse.h		\
	sparse.inl		\
	scalar-matrix.h		\
	submatrix.h		\
	inverse.h		\
	transpose.h		\
	sum.h			\
	diagonal.h		\
	moore-penrose.h		\
	dense.h  		\
	dense.inl		\
	direct-sum.h		\
	null-matrix.h		\
	matrix-blackbox.h	\
	dif.h			\
	companion.h		\
	frobenius.h		\
	zero-one.h		\
	zero-one.inl		\
	zo-gf2.h     zo-gf2.inl \
	triplesbb.h		\
	apply.h			\
	submatrix-traits.h	\
	random-matrix-traits.h	\
	scompose.h		\
	blackbox-interface.h	\
	random-matrix.h		\
	lambda-sparse.h		\
	blas-blackbox.h		\
	polynomial.h		\
	block-hankel.h		\
	block-toeplitz.h	\
	block-hankel-inverse.h	\
	jit-matrix.h		\
	rational-matrix-factory.h

NTL_HDRS =			\
	toeplitz.h ntl-hankel.h ntl-sylvester.h \
	toeplitz.inl ntl-hankel.inl ntl-sylvester.inl

GIVARO_HDRS = 

LIDIA_HDRS = 

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if LINBOX_HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

pkgincludesub_HEADERS =		\
	$(BASIC_HDRS)		\
	$(USE_NTL_HDRS)		\
	$(USE_GIVARO_HDRS)	\
	$(USE_LIDIA_HDRS)

EXTRA_DIST = 			\
	$(NTL_HDRS)		\
	$(GIVARO_HDRS)		\
	$(LIDIA_HDRS)		\
	blackbox.doxy		

