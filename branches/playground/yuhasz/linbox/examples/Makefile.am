INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS = $(GMP_CFLAGS) $(NTL_CFLAGS)
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS)

SUBDIRS=blackbox fields integers 
# not ready SUBDIRS+= gen_superlu

EXTRA_PROGRAMS=mmZZ_pE dot-product

mmZZ_pE_SOURCES = \
	mmZZ_pE.C

dot_product_SOURCES = \
	dot-product.C

LINBOX=@prefix@

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=@bindir@

# for compilation of new examples
%:%.C 
	$(CXX) $(CPPFLAGS) $*.C -o $*  $(shell $(LINBOX_BIN)/linbox-config --cflags --libs)
