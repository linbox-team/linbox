INCLUDES=-I$(top_srcdir)

SUBDIRS=data 

noinst_PROGRAMS = example ex_mat_0.C

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

example_SOURCES = \
	example.C

example_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

ex_mat_0_SOURCES = ex_mat_0.C 
ex_mat_0_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la
