SUBDIRS=data

noinst_PROGRAMS = example

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(GMP_LIBS) $(NTL_LIBS) $(top_builddir)/linbox/liblinbox.la

example_SOURCES = \
	example.C

example_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la