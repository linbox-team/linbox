# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/matrix

BASIC_HDRS =			\
	archetype.h		\
	matrix-category.h	\
	sparse.h		\
	sparse.inl		\
	transpose-matrix.h		\
	dense.h  		\
	dense.inl		\
	dense-submatrix.h	\
	dense-submatrix.inl	\
	dense-rows-matrix.h	\
	matrix-domain.h		\
	matrix-domain.inl	\
	blas-matrix.h		\
	factorized-matrix.h	\
	factorized-matrix.inl	

NTL_HDRS =
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
	matrix.doxy

