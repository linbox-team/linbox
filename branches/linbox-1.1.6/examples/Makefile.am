
#INCLUDES=-I$(top_srcdir)
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/linbox

#if LinBoxSrcOnly
#AM_CPPFLAGS=-g -DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
#LINBOX_LIBS=
#else
#AM_CPPFLAGS=-g $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
#LINBOX_LIBS= $(top_builddir)/linbox/liblinbox.la
#endif

AM_CXXFLAGS = -g -DLinBoxSrcOnly -Wall -O2
AM_CPPFLAGS+= $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
LDADD = $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS)

#EX_LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS)
#LDADD = $(LINBOX_LIBS) $(EX_LDADD)

SUBDIRS=fields

EXAMPLES=rank det minpoly valence solve dotproduct echelon sparseelimdet sparseelimrank checksolve doubledet
GIVARONTL_EXAMPLES=smith charpoly graph-charpoly
if LINBOX_HAVE_NTL
if LINBOX_HAVE_GIVARO
USE_GIVARONTL_EXAMPLES = $(GIVARONTL_EXAMPLES)
endif
endif


EXTRA_PROGRAMS= $(EXAMPLES) $(USE_GIVARONTL_EXAMPLES)

examples: ${EXTRA_PROGRAMS}

dotproduct_SOURCES = \
	dot-product.C

doubledet_SOURCES= doubledet.C
#linbox_program_SOURCES = \
#	linbox-program.C


checksolve_SOURCES = checksolve.C
charpoly_SOURCES = charpoly.C Matio.h
graph_charpoly_SOURCES = graph-charpoly.C Matio.h
det_SOURCES = det.C
rank_SOURCES = rank.C
smith_SOURCES = smith.C
minpoly_SOURCES = minpoly.C
valence_SOURCES = valence.C
solve_SOURCES = solve.C
echelon_SOURCES=echelon.C
sparseelimdet_SOURCES=sparseelimdet.C
sparseelimrank_SOURCES=sparseelimrank.C
LINBOX=@prefix@

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=@bindir@

# for compilation of new examples
%:%.C
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDADD) $(top_srcdir)/linbox/.libs/liblinbox.a
