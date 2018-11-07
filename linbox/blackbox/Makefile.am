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



pkgincludesubdir=$(pkgincludedir)/blackbox

BASIC_HDRS =			\
	apply.h                   \
	archetype.h               \
	bb.h                      \
	blackbox.h                \
	blackbox-interface.h      \
	blockbb.h                 \
	block-hankel.h            \
	block-hankel-inverse.h    \
	block-compose.h           \
	block-toeplitz.h          \
	butterfly.h               \
	butterfly.inl             \
	companion.h               \
	compose.h                 \
	csf.h                     \
	csf.inl                   \
	diagonal-gf2.h            \
	diagonal.h                \
	dif.h                     \
	direct-sum.h              \
	factory.h                 \
	fflas-csr.h               \
	fibb.h			          \
	fibb-product.h            \
	frobenius.h               \
	hilbert.h                 \
	inverse.h                 \
	jit-matrix.h              \
	lambda-sparse.h           \
	matrix-blackbox.h         \
	moore-penrose.h           \
	null-matrix.h             \
	pascal.h		          \
	permutation.h             \
	polynomial.h              \
	quad-matrix.h             \
	random-matrix.h           \
	random-matrix-traits.h    \
	rational-matrix-factory.h \
	scalar-matrix.h           \
	scompose.h                \
	squarize.h                \
	submatrix.h               \
	submatrix-traits.h        \
	sum.h                     \
	toeplitz.h                \
	toeplitz.inl              \
	transpose.h               \
	triangular-fibb.h         \
	zero-one.h                \
	zero-one.inl              \
	zo-gf2.h                  \
	zo-gf2.inl                \
	zo.h                      \
	zoi.inl                   \
	zo.inl

NTL_HDRS =			\
	ntl-hankel.h            \
	ntl-hankel.inl          \
	ntl-sylvester.h         \
	ntl-sylvester.inl


if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

pkgincludesub_HEADERS =		\
	$(BASIC_HDRS)		\
	$(USE_NTL_HDRS)		

EXTRA_DIST = 			\
	$(NTL_HDRS)		\
	blackbox.doxy

