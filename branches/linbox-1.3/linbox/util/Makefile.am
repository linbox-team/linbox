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



#we now need to include givaro headers for timer ? (oupa ? plus moche que ça tu meurs)
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/linbox $(FFLAS_FFPACK_CFLAGS) $(GIVARO_CFLAGS)

LIBS=$(GIVARO_LIBS)

AM_CXXFLAGS=@DEFAULT_CFLAGS@

SUBDIRS=formats

pkgincludesubdir=$(pkgincludedir)/util

noinst_LTLIBRARIES=libutil.la


libutil_la_SOURCES=	\
	error.C		    \
	debug.C

libutil_la_LDFLAGS=$(GIVARO_LIBS)  $(LDFLAGS)

#libutil_la_LIBADD=	\
#	gmp++/libgmpxx.la

pkgincludesub_HEADERS=      \
	debug.h				    \
	timer.h				    \
	error.h				    \
	commentator.h 			\
	commentator.inl 			\
	contracts.h 			\
	prime-stream.h			\
	iml_wrapper.h           \
	field-axpy.h			\
	matrix-stream.h			\
	matrix-stream.inl		\
	mpicpp.h			    \
	mpicpp.inl			\
	write-mm.h

EXTRA_DIST = util.doxy
