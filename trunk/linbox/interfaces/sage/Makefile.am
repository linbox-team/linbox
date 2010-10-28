# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
if LINBOX_HAVE_SAGE

INCLUDES=-I$(top_srcdir) -I. -I../../linbox

AM_CXXFLAGS = -O2 -DDISABLE_COMMENTATOR $(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(GIVARO_CFLAGS)

lib_LTLIBRARIES = liblinboxsage.la

liblinboxsage_la_SOURCES = linbox-sage.C

pkginclude_HEADERS = linbox-sage.h

#liblinboxsage_la_LIBADD = -llinbox $(BLAS_LIBS)

liblinboxsage_la_LDFLAGS =  $(GIVARO_LIBS) $(GMP_LIBS) $(NTL_LIBS) $(BLAS_LIBS) $(MAPLE_LIBS)  -version-info 0:0:0 #-Wl,-zmuldefs
endif