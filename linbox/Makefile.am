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

AM_CPPFLAGS=-I$(top_srcdir)/linbox
AM_CXXFLAGS=$(LINBOX_CXXFLAGS) $(FFLAS_FFPACK_CFLAGS) $(LINBOX_DEPS_CFLAGS)
LDADD=$(LINBOX_DEPS_LIBS) $(FFLAS_FFPACK_LIBS) $(LINBOX_LDFLAGS)

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
    $(LINBOX_DEPS_LIBS)

#  $(top_srcdir)/benchmarks/libbenchmarks.la \
#
