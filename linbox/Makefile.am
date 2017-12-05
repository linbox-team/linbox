# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
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
#/

AM_CPPFLAGS= -I$(top_srcdir)/linbox $(DEPS_CFLAGS) $(DEFAULT_CFLAGS)
LDADD = $(DEPS_LIBS) $(LDFLAGS)

#AM_CPPFLAGS = -I@includedir@ -I$(top_srcdir)/linbox -I$(top_srcdir) $(DEFCPPFLAGS) $(OPTFLAGS) $(FFLAS_FFPACK_CFLAGS) $(NTL_CFLAGS) $(MPFR_CFLAGS) $(FPLLL_CFLAGS) $(IML_CFLAGS) $(FLINT_CFLAGS) $(PARFLAGS)
#LDADD = $(FFLAS_FFPACK_LIBS) $(NTL_LIBS) $(MPFR_LIBS) $(FPLLL_LIBS) $(IML_LIBS) $(FLINT_LIBS) $(OCL_LIBS)

SUBDIRS=util element field randiter vector matrix blackbox algorithms solutions ring polynomial


BASIC_HDRS =          \
    linbox-config.h \
    config.h \
    config-blas.h     \
    integer.h         \
    linbox-tags.h

NTL_HDRS =

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

pkginclude_HEADERS =        \
    $(BASIC_HDRS)           \
    $(USE_NTL_HDRS)

EXTRA_DIST =            \
    $(NTL_HDRS)         \
    linbox.doxy         \
    archetypes.doxy

lib_LTLIBRARIES=liblinbox.la

liblinbox_la_SOURCES = dummy.C

liblinbox_la_LIBADD =           \
    util/libutil.la             \
    randiter/libranditer.la     \
    algorithms/libalgorithms.la \
    $(MPFR_LIBS) $(IML_LIBS) $(OCL_LIBS)

#Cygwin ?
liblinbox_la_LDFLAGS = $(DEPS_LIBS) $(LDFLAGS) -no-undefined


#  $(top_srcdir)/benchmarks/libbenchmarks.la \
#
