SUBDIRS=util element field randiter vector blackbox algorithms solutions

EXTRA_DIST=linbox.dxx objects.dxx

includedir=$(prefix)/include/linbox

include_HEADERS=			\
	debug.h				\
	error.h				\
	integer.h	integer.inl

lib_LTLIBRARIES=liblinbox.la

liblinbox_la_SOURCES=	\
	error.C		\
	debug.C

liblinbox_la_LIBADD= \
	util/libutil.la
