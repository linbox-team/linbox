INCLUDES=-I$(top_srcdir)

SUBDIRS=

noinst_PROGRAMS = first

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

first_SOURCES = \
	first.C

first_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la
