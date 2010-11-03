# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence

#INCLUDES=-I$(top_srcdir)
INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/linbox

#if LinBoxSrcOnly
#AM_CPPFLAGS=-g -DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) 
#LINBOX_LIBS= 
#else
#AM_CPPFLAGS=-g $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) 
#LINBOX_LIBS= $(top_builddir)/linbox/liblinbox.la 
#endif

#SAFER_FLAGS=--std=c++98 -g -Wall -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter
SAFER_FLAGS=-g -Wall

AM_CXXFLAGS = $(SAFER_FLAGS) -DLinBoxSrcOnly -DDEBUG -O0 
# AM_CXXFLAGS = $(SAFER_FLAGS) -DLinBoxSrcOnly -O2 
AM_CPPFLAGS= $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
#LDADD = $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS)  "-Wl,-enable-auto-import"
LDADD = $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS) 

#EX_LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) 
#LDADD = $(LINBOX_LIBS) $(EX_LDADD)

SUBDIRS=fields 

EXAMPLES=rank det minpoly valence solve dot-product echelon sparseelimdet sparseelimrank checksolve doubledet smithvalence
# EXAMPLES+=nulp yabla
GIVARONTL_EXAMPLES=smith charpoly graph-charpoly
if LINBOX_HAVE_NTL
if LINBOX_HAVE_GIVARO
USE_GIVARONTL_EXAMPLES = $(GIVARONTL_EXAMPLES)
endif
endif


EXTRA_PROGRAMS= $(EXAMPLES) $(USE_GIVARONTL_EXAMPLES) 

examples: ${EXTRA_PROGRAMS}

dot_product_SOURCES = \
	dot-product.C

doubledet_SOURCES= doubledet.C
#linbox_program_SOURCES = \
#	linbox-program.C


checksolve_SOURCES = checksolve.C
charpoly_SOURCES = charpoly.C
graph_charpoly_SOURCES = graph-charpoly.C
det_SOURCES = det.C
rank_SOURCES = rank.C
smith_SOURCES = smith.C
minpoly_SOURCES = minpoly.C
valence_SOURCES = valence.C
solve_SOURCES = solve.C
echelon_SOURCES=echelon.C
smithvalence_SOURCES=smithvalence.C
sparseelimdet_SOURCES=sparseelimdet.C
sparseelimrank_SOURCES=sparseelimrank.C
LINBOX=@prefix@

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=@bindir@


# for compilation of new examples
%:%.C 
	$(CXX) $(SAFER_FLAGS) $(CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(top_srcdir)/linbox/.libs/liblinbox.a

%:%.cpp
	$(CXX) $(SAFER_FLAGS) $(CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.cpp -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(top_srcdir)/linbox/.libs/liblinbox.a
