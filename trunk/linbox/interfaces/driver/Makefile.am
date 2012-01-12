# Copyright (c) 2010 the LinBox group
#========LICENCE========
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
# License along with this library; if not, write to the Free Software Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========

if LINBOX_COMPILE_DRIVERS
pkgincludesubdir=$(includedir)/interfaces/driver

AM_CPPFLAGS=-I$(top_srcdir) -I. -I../../linbox

AM_CXXFLAGS = @DEFAULT_CFLAGS@ -DDISABLE_COMMENTATOR $(GMP_CFLAGS) $(FFLAS_FFPACK_CFLAGS) $(BLAS_CFLAGS) $(GIVARO_CFLAGS)

#ifdef HAVE_NTL
AM_CXXFLAGS +=$(NTL_CFLAGS)
#endif

#LDADD= $(GMP_LIBS) $(NTL_LIBS) $(BLAS_LIBS) $(MAPLE_LIBS) $(GIVARO_LIBS) "../../linbox/liblinbox.la -z muldefs"

lib_LTLIBRARIES=liblbdriver.la


liblbdriver_la_SOURCES= lb-element.C 	\
			lb-domain.C     \
			lb-vector.C     \
			lb-polynomial.C \
			lb-det.C        \
			lb-rank.C       \
			lb-garbage.C    \
			lb-minpoly.C    \
			lb-blackbox.C   \
			lb-charpoly.C
# \
			#  lb-solve.C

liblbdriver_la_LDFLAGS=$(GIVARO_LIBS) $(GMP_LIBS) $(NTL_LIBS) $(BLAS_LIBS) $(LDFLAGS) $(top_srcdir)/linbox/liblinbox.la -Wl,-zmuldefs



pkginclude_HEADERS=\
		lb-driver.h              \
		lb-blackbox-abstract.h   \
		lb-blackbox-type.h       \
		lb-domain-function.h     \
		lb-element-collection.h  \
		lb-polynomial.h          \
		lb-vector-data.h         \
		lb-blackbox-collection.h \
		lb-charpoly.h            \
		lb-domain-functor.h      \
		lb-element-data.h        \
		lb-rank.h                \
		lb-vector-function.h     \
		lb-blackbox-data.h       \
		lb-det.h                 \
		lb-domain.h              \
		lb-element.h             \
		lb-solve.h               \
		lb-vector-functor.h      \
		lb-blackbox-function.h   \
		lb-domain-abstract.h     \
		lb-domain-type.h         \
		lb-garbage.h             \
		lb-utilities.h           \
		lb-vector.h              \
		lb-blackbox-functor.h    \
		lb-domain-collection.h   \
		lb-interface.h           \
		lb-vector-abstract.h     \
		lb-vector-type.h         \
		lb-blackbox.h            \
		lb-domain-data.h         \
		lb-element-abstract.h    \
		lb-minpoly.h             \
		lb-vector-collection.h
endif
