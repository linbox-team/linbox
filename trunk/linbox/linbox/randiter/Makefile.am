includedir=$(prefix)/include/linbox/randiter

EXTRA_DIST=		\
	archetype.dxx	\
	random.dxx

BASIC_HDRS =				\
	abstract.h			\
	archetype.h			\
	envelope.h			\
	gmp-rational.h			\
	param-fuzzy.h			\
	modular.h			\
	unparametric.h			\
	nonzero.h

NTL_HDRS =

GIVARO_HDRS =

LIDIA_HDRS =				\
	lidia.h				\
	lidia-gfq.h

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

EXTRA_HEADERS =				\
	$(NTL_HDRS)			\
	$(GIVARO_HDRS)			\
	$(LIDIA_HDRS)
