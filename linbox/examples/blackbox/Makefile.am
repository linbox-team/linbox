INCLUDES=-I$(top_srcdir)

noinst_PROGRAMS = example ex-mat0

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

example_SOURCES = \
	example.C

example_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

ex_mat0_SOURCES = ex-mat0.C
ex_mat0_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

# general 

%:%.o
	$(CXXLINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)



