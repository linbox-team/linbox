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


pkgincludesubdir=$(pkgincludedir)/matrix

SUBDIRS= SparseMatrix

BASIC_HDRS =              \
	archetype.h               \
	matrix-category.h         \
	matrix-traits.h           \
	sparse.h                  \
	transpose-matrix.h        \
	dense-rows-matrix.h       \
	matrix-domain.h           \
	matrix-domain.inl         \
	matrix-domain-gf2.h       \
	blas-matrix.h             \
	blas-matrix.inl           \
	blas-matrix-domain.h      \
	blas-matrix-domain.inl    \
	sparse-matrix.h           \
	sparse-formats.h          \
	blas-submatrix.inl        \
	blas-matrix-multimod.h    \
	factorized-matrix.h       \
	factorized-matrix.inl     \
	blas-triangularmatrix.inl \
	permutation-matrix.h      \
	permutation-matrix.inl    \
	random-matrix.h           \
	random-matrix.inl

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

pkgincludesub_HEADERS =     \
    $(BASIC_HDRS)           \
    $(USE_NTL_HDRS)         \
    $(USE_GIVARO_HDRS)      \
    $(USE_LIDIA_HDRS)

EXTRA_DIST =            \
    $(NTL_HDRS)         \
    $(GIVARO_HDRS)      \
    $(LIDIA_HDRS)       \
    matrix.doxy

