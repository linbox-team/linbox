INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(MAPLE_CFLAGS) $(GIVARO_CFLAGS)

EXTRA_DIST = lbmaple.C lbmaple.mpl.head lbmaple.mpl.tail Makefile.in.1 Makefile.in.2 demonstration.mws

lib_LTLIBRARIES=liblbmaple.la liblbmapleffpack.la

liblbmaple_la_SOURCES=lbmaple.C
liblbmaple_la_LDFLAGS=$(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) $(MAPLE_LIBS)

liblbmapleffpack_la_SOURCES=lbmaple-ffpack.C
liblbmapleffpack_la_LDFLAGS=$(top_builddir)/linbox/liblinbox.la $(ATLAS_LIBS) $(GMP_LIBS) $(NTL_LIBS) $(MAPLE_LIBS) $(GIVARO_LIBS)

install-data-local:
	$(INSTALL_DATA) lbmaple.mpl $(libdir)/lbmaple.mpl
