INCLUDES=-I$(top_srcdir)

SUBDIRS=util element field randiter vector blackbox algorithms solutions

EXTRA_DIST=linbox.dxx objects.dxx

includedir=$(prefix)/include/linbox

include_HEADERS=			\
	integer.h	integer.inl	\
	$(NTL_HEADERS_BASE)		\
	$(GIVARO_HEADERS_BASE)

EXTRA_HEADERS=	ntl.h givaro.h

lib_LTLIBRARIES=liblinbox.la

liblinbox_la_SOURCES=dummy.C

liblinbox_la_LIBADD= \
	util/libutil.la
