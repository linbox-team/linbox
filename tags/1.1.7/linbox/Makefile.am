# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir)

SUBDIRS=util element field randiter vector matrix blackbox switch algorithms solutions fflas ffpack ring

BASIC_HDRS =				\
	integer.h	\
	config-blas.h	linbox-config.h

NTL_HDRS =

GIVARO_HDRS =

LIDIA_HDRS =

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if LINBOX_HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

pkginclude_HEADERS =			\
	$(BASIC_HDRS)			\
	$(USE_NTL_HDRS)			\
	$(USE_GIVARO_HDRS)		\
	$(USE_LIDIA_HDRS)

EXTRA_DIST =				\
	$(NTL_HDRS)			\
	$(GIVARO_HDRS)			\
	$(LIDIA_HDRS)			\
	linbox.doxy			\
	archetypes.doxy

lib_LTLIBRARIES=liblinbox.la

liblinbox_la_SOURCES = dummy.C

liblinbox_la_LIBADD = \
	util/libutil.la		\
	randiter/libranditer.la \
	algorithms/libalgorithms.la
