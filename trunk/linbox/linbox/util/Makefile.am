# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
#we now need to include givaro headers for timer ? (oupa ? plus moche que ça tu meurs)
INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/linbox $(GIVARO_CFLAGS)
LIBS=-lgivaro



AM_CXXFLAGS=@DEFAULT_CFLAGS@

SUBDIRS=formats

pkgincludesubdir=$(pkgincludedir)/util

noinst_LTLIBRARIES=libutil.la


libutil_la_SOURCES=	\
	error.C		    \
	commentator.C	\
	debug.C

#libutil_la_LIBADD=	\
#	gmp++/libgmpxx.la

pkgincludesub_HEADERS=      \
	debug.h				    \
	timer.h				    \
	error.h				    \
	commentator.h 			\
	contracts.h 			\
	prime-stream.h			\
	field-axpy.h			\
	matrix-stream.h			\
	matrix-stream.inl		\
	mpicpp.h			    \
	mpicpp.inl

EXTRA_DIST = util.doxy
