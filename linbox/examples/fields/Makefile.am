INCLUDES=-I$(top_srcdir)

noinst_PROGRAMS = ex-fields.C  ex-fields-0.C

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 


ex-fields-0_SOURCES = ex-fields-0.C 
ex-fields-0_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

ex-fields_SOURCES = ex-fields.C 
ex-fields_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la
