INCLUDES=-I$(top_srcdir)

noinst_PROGRAMS = ex-fields0 ex-fields

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 


ex_fields0_SOURCES = ex-fields0.C 
ex_fields0_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la


