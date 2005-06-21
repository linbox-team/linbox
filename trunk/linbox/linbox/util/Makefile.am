INCLUDES=-I$(top_srcdir)

SUBDIRS=gmp++ 

pkgincludesubdir=$(pkgincludedir)/util

noinst_LTLIBRARIES=libutil.la

libutil_la_SOURCES=	\
	commentator.C	\
	timer.C		\
	error.C		\
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
