INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS = $(GMP_CFLAGS) $(NTL_CFLAGS)
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS)

#SUBDIRS=

EXTRA_PROGRAMS=genrank

genrank_SOURCES = \
	genrank.C

