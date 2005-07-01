
INCLUDES=-I$(top_srcdir)

SUBDIRS=gmp++ formats

pkgincludesubdir=$(pkgincludedir)/util

noinst_LTLIBRARIES=libutil.la


libutil_la_SOURCES=	\
	error.C		\
	commentator.C	\
	debug.C

libutil_la_LIBADD=	\
	gmp++/libgmpxx.la

pkgincludesub_HEADERS=\
	debug.h				\
	error.h				\
	commentator.h 			\
	timer.h				\
	prime-stream.h			\
	field-axpy.h
EXTRA_DIST = util.doxy
