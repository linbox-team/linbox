INCLUDES=-I$(top_srcdir)

EXTRA_PROGRAMS = smithbig smith example ex-mat0 load-det load-minpoly load-symmetrize random-rank 2localTref smithhybrid

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) -lgivaro -lclapack -lcblas -latlas

example_SOURCES  = example.C
smith_SOURCES  = smith.C
2localTref_SOURCES  = 2localTref.C
smithbig_SOURCES  = smithbig.C
load_det_SOURCES = load-det.C
load_minpoly_SOURCES = load-minpoly.C
load_symmetrize_SOURCES = load-symmetrize.C
random_rank_SOURCES = random-rank.C
ex_mat0_SOURCES  = ex-mat0.C
smithhybrid_SOURCES  = smithhybrid.C

# general 

%:%.o
	$(CXXLINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)



