INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) @MAPLE_INCLUDES@ 
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 
EXTRA_DIST = lbmaple.C lbmaple.mpl.head lbmaple.mpl.tail Makefile.blank maplebb.h demonstration.mws

mapledir=$(top_builddir)/linbox/util/maple
maple_LTLIBRARIES=liblbmaple.la
libmaple_la_SOURCES=lbmaple.C
libmaple_la_INCLUDES=maplebb.h
libmaple_la_LDFLAGS=@MAPLE_LIBS@ -module

