includedir=$(prefix)/include/linbox/field

BASIC_HDRS =				\
	field-interface.h		\
	abstract.h			\
	archetype.h			\
	envelope.h			\
	unparametric.h			\
	param-fuzzy.h			\
	modular.h	modular.inl	\
	gmp-rational.h			\
	givaro-zpz.h 			\
	FieldBLAS.h	FieldBLAS.C	\
	vector-domain.h	vector-domain.C

NTL_HDRS =				\
	ntl.h				\
	ntl-RR.h			\
	ntl-ZZ_p.h			\
	ntl-ZZ_pE.h			\
	ntl-zz_p.h

GIVARO_HDRS =				\
	givaro.h			\
	givaro-zpz.h			\
	givaro-gfq.h

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

EXTRA_DIST=				\
	abstract.dxx			\
	fields.dxx			\
	parametric.dxx			\
	unparametric.dxx		\
	wrappers.dxx
