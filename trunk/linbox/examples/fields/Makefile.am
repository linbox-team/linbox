INCLUDES=-I$(top_srcdir)

noinst_PROGRAMS = ex-fields0 ex-fields \
			ex-fields-wrapper  ex-fields-archetype

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 


ex_fields_SOURCES = ex-fields.C 
ex_fields_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

ex_fields0_SOURCES = ex-fields0.C 
ex_fields0_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

ex_fields_wrapper_SOURCES = ex-fields-wrapper.C
ex_fields_wrapper_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la

ex_fields_archetype_SOURCES = ex-fields-archetype.C
ex_fields_archetype_DEPENDENCIES = $(top_builddir)/linbox/liblinbox.la


