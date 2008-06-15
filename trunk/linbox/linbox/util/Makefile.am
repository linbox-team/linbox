INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/linbox

AM_CFLAGS=-Wall -O2

SUBDIRS=formats

pkgincludesubdir=$(pkgincludedir)/util

noinst_LTLIBRARIES=libutil.la


libutil_la_SOURCES=	\
	timer.C		\
	error.C		\
	commentator.C	\
	debug.C

#libutil_la_LIBADD=	\
#	gmp++/libgmpxx.la

pkgincludesub_HEADERS=\
	debug.h				\
	error.h				\
	commentator.h 			\
	timer.h				\
	prime-stream.h			\
	field-axpy.h			\
	matrix-stream.h			\
	matrix-stream.inl		\
	mpicpp.h			\
	mpicpp.inl

EXTRA_DIST = util.doxy
