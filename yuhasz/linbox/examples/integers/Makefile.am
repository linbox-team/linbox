INCLUDES=-I$(top_srcdir)

EXTRA_PROGRAMS = ex-integers

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

ex_integers_SOURCES = ex-integers.C


