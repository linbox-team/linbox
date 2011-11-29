if LINBOX_HAVE_SAGE

AM_CPPFLAGS=-I$(top_srcdir) -I. -I../../linbox

AM_CXXFLAGS = -O2 -DDISABLE_COMMENTATOR $(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(GIVARO_CFLAGS)

lib_LTLIBRARIES = liblinboxsage.la

liblinboxsage_la_SOURCES = linbox-sage.C

pkginclude_HEADERS = linbox-sage.h

liblinboxsage_la_LIBADD = $(top_builddir)/linbox/liblinbox.la

liblinboxsage_la_LDFLAGS =  $(GIVARO_LIBS) $(GMP_LIBS) $(NTL_LIBS) $(BLAS_LIBS) $(MAPLE_LIBS)  -version-info 0:0:0 #-Wl,-zmuldefs
endif
