# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
if LINBOX_HAVE_MAPLE
pkgincludesubdir=$(includedir)/interfaces/maple/new

INCLUDES=-I$(top_srcdir) -I. -I../../linbox

AM_CXXFLAGS = -O2 -DDISABLE_COMMENTATOR -I../driver $(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(MAPLE_CFLAGS) $(GIVARO_CFLAGS)

lib_LTLIBRARIES=liblbmaple.la 

liblbmaple_la_SOURCES= lb-maple.C

liblbmaple_la_LDFLAGS=$(GIVARO_LIBS) $(GMP_LIBS) $(NTL_LIBS) $(BLAS_LIBS) $(MAPLE_LIBS) $(top_srcdir)/interfaces/driver/liblbdriver.la -Wl,-zmuldefs

pkginclude_HEADERS =lb-maple.h lb-maple-utilities.h


install-data-local: 
	sh lb-maple-path.sh "$(libdir)/liblbmaple.so"
	$(INSTALL_DATA) lb-maple.mpl $(libdir)/lb-maple.mpl

install-data-hook:
	cd $(libdir) && \
	$(MAPLE_HOME)/bin/maple lb-maple.mpl 2>&1 > /dev/null
endif