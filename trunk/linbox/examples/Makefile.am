INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=-g $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(ATLAS_CFLAGS) 
LDADD = $(top_builddir)/linbox/liblinbox.la $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(ATLAS_LIBS) 

SUBDIRS=blackbox fields integers 
# not ready SUBDIRS+= gen_superlu

EXTRA_PROGRAMS=mmZZ_pE dot-product dense-charpoly

mmZZ_pE_SOURCES = \
	mmZZ_pE.C

dot_product_SOURCES = \
	dot-product.C

dense_charpoly_SOURCES = \
	dense-charpoly.C

LINBOX=@prefix@

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=@bindir@

# for compilation of new examples
%:%.C 
	$(CXX) $(CPPFLAGS) $*.C -o $*  $(shell $(LINBOX_BIN)/linbox-config --cflags --libs)
