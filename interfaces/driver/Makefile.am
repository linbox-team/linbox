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
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========

if LINBOX_COMPILE_DRIVERS
pkgincludesubdir=$(includedir)/interfaces/driver

AM_CPPFLAGS=-I$(top_srcdir) -I. -I../../linbox -DDISABLE_COMMENTATOR
AM_CXXFLAGS=$(LINBOX_CXXFLAGS) $(FFLAS_FFPACK_CFLAGS) $(LINBOX_DEPS_CFLAGS)

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

liblbdriver_la_LDFLAGS=$(top_srcdir)/linbox/liblinbox.la $(LINBOX_DEPS_LIBS) $(FFLAS_FFPACK_LIBS) $(LINBOX_LDFLAGS) -Wl,-zmuldefs


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
