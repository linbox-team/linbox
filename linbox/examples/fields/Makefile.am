INCLUDES=-I$(top_srcdir)

noinst_PROGRAMS = ex_fields.C  ex_fields_0.C

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 


ex_fields_0_SOURCES = ex_fields_0.C 
ex_fields_0_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

ex_fields_SOURCES = ex_fields.C 
ex_fields_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la
