# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=-g $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(ATLAS_CFLAGS)
LDADD = $(top_builddir)/linbox/liblinbox.la $(NTL_LIBS) $(GIVARO_LIBS)  $(GMP_LIBS) $(LIDIA_LIBS) $(ATLAS_LIBS) $(LDFLAGS)

EXTRA_PROGRAMS=t-rdisolve

t_rdisolve_SOURCES = \
	t-rdisolve.C

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=$(bindir)

# for compilation of new examples
%:%.C
	$(CXX) $(CPPFLAGS) $*.C -o $*  $(shell $(LINBOX_BIN)/linbox-config --cflags --libs) $(LDADD)
