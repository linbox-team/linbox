INCLUDES=-I$(top_srcdir)

SUBDIRS=integer

includedir=$(prefix)/include/linbox/util

noinst_LTLIBRARIES=libutil.la

libutil_la_SOURCES=	\
	commentator.C	\
	timer.C

libutil_la_LIBADD=	\
	integer/libinteger.la

include_HEADERS=\
	commentator.h	commentator.inl	\
	timer.h				\
	field-axpy.h