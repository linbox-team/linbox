INCLUDES=-I$(top_srcdir)

SUBDIRS=

noinst_PROGRAMS =

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 


# general 

%:%.o
	$(CXXLINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)


