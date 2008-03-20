INCLUDES=-I$(top_srcdir)

EXTRA_PROGRAMS = smith example ex-mat0 load-det load-minpoly load-symmetrize random-rank 
#zeroone-rank

AM_CPPFLAGS=-DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(GMP_LIBS) $(NTL_LIBS) -lgivaro -lclapack -lcblas -latlas 

example_SOURCES  = example.C
smith_SOURCES  = smith.C
load_det_SOURCES = load-det.C
load_minpoly_SOURCES = load-minpoly.C
load_symmetrize_SOURCES = load-symmetrize.C
random_rank_SOURCES = random-rank.C
#zeroone_rank_SOURCES = zeroone-rank.C
ex_mat0_SOURCES  = ex-mat0.C

# general 

%:%.o
	$(CXXLINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)


