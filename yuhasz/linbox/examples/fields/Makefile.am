INCLUDES=-I$(top_srcdir)

EXTRA_PROGRAMS =		\
	ex-fields0		\
	ex-fields 		\
	ex-fields-wrapper	\
	ex-fields-archetype	\
	modular-int

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) 

ex_fields_SOURCES           = ex-fields.C 
ex_fields0_SOURCES          = ex-fields0.C 
ex_fields_wrapper_SOURCES   = ex-fields-wrapper.C
ex_fields_archetype_SOURCES = ex-fields-archetype.C
modular_int_SOURCES = modular-int.C
