# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/field

BASIC_HDRS =			\
	field-interface.h	\
	field-traits.h		\
	rebind.h		\
	abstract.h		\
	archetype.h		\
	envelope.h		\
	unparametric.h		\
	param-fuzzy.h		\
	gmp-rational.h		\
	modular.h		\
	modular.inl		\
	modular-int.h		\
	modular-int32.h		\
	modular-short.h		\
	modular-byte.h		\
	modular-balanced-double.h	\
	modular-balanced-float.h	\
	modular-balanced-int.h		\
	modular-balanced-int32.h	\
	modular-double.h	\
	modular-float.h		\
	gf2.h			\
	gf2.inl			\
	hom.h			\
	map.h			\
	PIR-modular-int.h	\
	PIR-modular-int32.h	\
	local2_32.h		\
	PID-integer.h		\
	PID-double.h		\
	block-ring.h		\
	multimod-field.h

NTL_HDRS =				\
	ntl.h				\
	ntl-RR.h			\
	ntl-ZZ_p.h			\
	ntl-ZZ_pE.h			\
	ntl-ZZ_pX.h			\
	ntl-lzz_p.h			\
	ntl-ZZ.h			\
	ntl-lzz_pE.h			\
	ntl-lzz_pX.h			\
	ntl-GF2E.h			\
	ntl-pid-lzz_p.h			\
	PIR-ntl-ZZ_p.h			

GIVARO_HDRS =				\
	givaro.h			\
	givaro-zpz.h			\
	givaro-montg.h			\
	givaro-zpz.inl			\
	givaro-gfq.h			\
	givaro-extension.h		\
	givaro-field.h		\
	givaro-rational.h

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
	field.doxy
