pkgincludesubdir=$(pkgincludedir)/field

BASIC_HDRS =			\
	field-interface.h	\
	field-traits.h		\
	abstract.h		\
	archetype.h		\
	envelope.h		\
	unparametric.h		\
	param-fuzzy.h		\
	gmp-rational.h		\
	givaro-zpz.h 		\
	givaro-montg.h 		\
	givaro-zpz.inl 		\
	modular.h		\
	modular.inl		\
	modular-int.h		\
	modular-short.h		\
	modular-byte.h		\
	modular-balance-int.h	\
	modular-double.h	\
	gf2.h			\
	gf2.inl			\
	hom.h			\
	image-field.h		\
	PIR-ntl-ZZ_p.h		\
	PIR-modular-int.h	\
	ntl-pid-zz_p.h		\
	local2_32.h

NTL_HDRS =				\
	ntl.h				\
	ntl-RR.h			\
	ntl-ZZ_p.h			\
	ntl-ZZ_pE.h			\
	ntl-zz_p.h			\
	ntl-ZZ.h

GIVARO_HDRS =				\
	givaro.h			\
	givaro-zpz.h			\
	givaro-montg.h			\
	givaro-zpz.inl			\
	givaro-gfq.h

LIDIA_HDRS =				\
	lidia.h				\
	lidia-gfq.h

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if LINBOX_HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

pkgincludesub_HEADERS =			\
	$(BASIC_HDRS)			\
	$(USE_NTL_HDRS)			\
	$(USE_GIVARO_HDRS)		\
	$(USE_LIDIA_HDRS)	

EXTRA_DIST =				\
	$(NTL_HDRS)			\
	$(GIVARO_HDRS)			\
	$(LIDIA_HDRS)			\
	abstract.dxx			\
	field.dxx			\
	parametric.dxx			\
	unparametric.dxx		\
	wrappers.dxx
