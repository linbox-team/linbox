INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS = $(GMP_CFLAGS) $(NTL_CFLAGS)
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS)

SUBDIRS=blackbox fields integers

EXTRA_PROGRAMS=mmZZ_pE

mmZZ_pE_SOURCES = \
	mmZZ_pE.C
