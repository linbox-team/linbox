INCLUDES=-I$(top_srcdir)

EXTRA_PROGRAMS = example ex-mat0 load-det

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

example_SOURCES  = example.C
load_det_SOURCES = load-det.C
ex_mat0_SOURCES  = ex-mat0.C

# general 

%:%.o
	$(CXXLINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)



