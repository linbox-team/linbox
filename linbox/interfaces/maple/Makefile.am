INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) $(MAPLE_CFLAGS)

EXTRA_DIST = lbmaple.C lbmaple.mpl.head lbmaple.mpl.tail Makefile.in.1 Makefile.in.2 demonstration.mws

lib_LTLIBRARIES=liblbmaple.la

liblbmaple_la_SOURCES=lbmaple.C
liblbmaple_la_LDFLAGS=$(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) $(MAPLE_LIBS)

install-data-local:
	$(INSTALL_DATA) lbmaple.mpl $(libdir)/lbmaple.mpl
