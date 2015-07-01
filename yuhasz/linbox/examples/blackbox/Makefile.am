INCLUDES=-I$(top_srcdir)

EXTRA_PROGRAMS = example ex-mat0 load-det load-minpoly load-symmetrize random-rank

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) -lgivaro

example_SOURCES  = example.C
load_det_SOURCES = load-det.C
load_minpoly_SOURCES = load-minpoly.C
load_symmetrize_SOURCES = load-symmetrize.C
random_rank_SOURCES = random-rank.C
ex_mat0_SOURCES  = ex-mat0.C

# general 

%:%.o
	$(CXXLINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)



