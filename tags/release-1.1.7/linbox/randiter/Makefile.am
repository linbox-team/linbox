# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir)/linbox $(GMP_CFLAGS)
pkgincludesubdir=$(pkgincludedir)/randiter

AM_CXXFLAGS = -O2 -Wall -DDISABLE_COMMENTATOR $(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(GIVARO_CFLAGS)
 
BASIC_HDRS =				\
	abstract.h			\
	archetype.h			\
	envelope.h			\
	gmp-rational.h			\
	param-fuzzy.h			\
	modular.h			\
	modular-balanced.h		\
	unparametric.h			\
	gf2.h				\
	nonzero.h			\
	mersenne-twister.h		\
	ntl-ZZ.h			\
	random-prime.h			\
	gmp-random-prime.h		\
	random-fftprime.h		\
	multimod-randomprime.h

NTL_HDRS =

GIVARO_HDRS =

LIDIA_HDRS =				\
	lidia.h				\
	lidia-gfq.h

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if LINBOX_HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

pkgincludesub_HEADERS =			\
	$(BASIC_HDRS)			\
	$(USE_NTL_HDRS)			\
	$(USE_GIVARO_HDRS)		\
	$(USE_LIDIA_HDRS)

EXTRA_DIST =				\
	$(NTL_HDRS)			\
	$(GIVARO_HDRS)			\
	$(LIDIA_HDRS)			\
	randiter.doxy			\
	archetype.doxy			\
	randiter-wrappers.doxy

noinst_LTLIBRARIES = libranditer.la

libranditer_la_SOURCES =	\
	mersenne-twister.C

