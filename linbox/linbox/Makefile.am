INCLUDES=-I$(top_srcdir)

SUBDIRS=util element field randiter vector matrix blackbox switch algorithms solutions FFLAS

includedir=$(prefix)/include/linbox

BASIC_HDRS =				\
	integer.h	integer.inl	\
	rational.h

NTL_HDRS =				\
	ntl.h

GIVARO_HDRS =				\
	givaro.h

LIDIA_HDRS =				\
	lidia.h

if HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

if HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

include_HEADERS =			\
	$(BASIC_HDRS)			\
	$(USE_NTL_HDRS)			\
	$(USE_GIVARO_HDRS)		\
	$(USE_LIDIA_HDRS)

EXTRA_DIST =				\
	$(NTL_HDRS)			\
	$(GIVARO_HDRS)			\
	$(LIDIA_HDRS)			\
	linbox.dxx			\
	interfaces.dxx			\
	objects.dxx

lib_LTLIBRARIES=liblinbox.la

liblinbox_la_SOURCES = dummy.C

liblinbox_la_LIBADD = \
	util/libutil.la		\
	randiter/libranditer.la
