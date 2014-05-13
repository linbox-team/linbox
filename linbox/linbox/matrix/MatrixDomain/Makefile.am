# Copyright (c) 2013 the LinBox group
# written by BB
#  ========LICENCE========
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



pkgincludesubdir=$(pkgincludedir)/matrix/MatrixDomain


OCL_HDRS =                                  \
    opencl-domain.h                         \
    opencl-domain.inl                       \
    opencl-domain-util.inl                  \
    opencl-domain-memory.inl                \
    opencl-domain-factory.h

#  if LINBOX_HAVE_OCL // we include them anyway. (for make dist when OCL not available)
USE_OCL_HDRS = $(OCL_HDRS)
#  endif



pkgincludesub_HEADERS =           \
	matrix-domain.h           \
	matrix-domain.inl         \
	matrix-domain-gf2.h       \
	apply-domain.h            \
	blas-matrix-domain.h      \
	blas-matrix-domain.inl    \
	plain-domain.h            \
	apply-domain.h            \
	$(USE_OCL_HDRS)




