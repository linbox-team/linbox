INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) @MAPLE_INCLUDES@ 
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 
EXTRA_DIST = lbmaple.C lbmaple.mpl.head lbmaple.mpl.tail Makefile.in.1 Makefile.in.2 maplebb.h \
		demonstration.mws

if HAVE_MAPLE
mapledir=@MAPLE_BUILD_LOC@
maple_LTLIBRARIES=liblbmaple.la
liblbmaple_la_SOURCES=lbmaple.C
liblbmaple_la_INCLUDES=maplebb.h
liblbmaple_la_LDFLAGS=@MAPLE_LIBS@ -module
endif

