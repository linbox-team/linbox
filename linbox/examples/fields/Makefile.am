# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
INCLUDES=-I$(top_srcdir)

EXTRA_PROGRAMS =		\
	ex-fields 		\
	ex-fields-archetype	\
	modular-int

#AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
#LDADD = $(top_builddir)/linbox/liblinbox.a $(GMP_LIBS) $(NTL_LIBS) 
AM_CPPFLAGS= -DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(GMP_LIBS) $(NTL_LIBS) -lgivaro -lclapack -lcblas -latlas

ex_fields_SOURCES           = ex-fields.C 
ex_fields_archetype_SOURCES = ex-fields-archetype.C
modular_int_SOURCES = modular-int.C
