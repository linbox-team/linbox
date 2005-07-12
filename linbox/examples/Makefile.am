
INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=-g -DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(ATLAS_CFLAGS) 
LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(ATLAS_LIBS) -lcblas -latlas 
#LDADD = $(top_builddir)/linbox/liblinbox.la $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(ATLAS_LIBS) 

SUBDIRS=fields 

EXTRA_PROGRAMS= rank det minpoly smith charpoly solve dot-product

examples: ${EXTRA_PROGRAMS}

dot_product_SOURCES = \
	dot-product.C

#linbox_program_SOURCES = \
#	linbox-program.C

charpoly_SOURCES = charpoly.C Matio.h
det_SOURCES = det.C
rank_SOURCES = rank.C
smith_SOURCES = smith.C
minpoly_SOURCES = minpoly.C
solve_SOURCES = solve.C

LINBOX=@prefix@

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=@bindir@

# for compilation of new examples
%:%.C 
	$(CXX) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $*  $(shell $(LINBOX_BIN)/linbox-config --cflags ) $(LDADD)
