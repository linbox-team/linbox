includedir=$(prefix)/include/linbox/field

include_HEADERS=			\
	field-interface.h		\
	abstract.h			\
	archetype.h			\
	envelope.h			\
	unparametric.h			\
	param-fuzzy.h			\
	modular.h	modular.inl	\
	gmp-rational.h			\
	givaro-zpz.h 			\
	FieldBLAS.h  FieldBLAS.C        \
	vector-domain.h	vector-domain.C	\
	$(NTL_HEADERS_FIELD)		\
	$(GIVARO_HEADERS_FIELD)		\
	$(LIDIA_HEADERS_FIELD)	

EXTRA_HEADERS= ntl.h ntl-RR.h ntl-ZZ_p.h ntl-ZZ_pE.h ntl-zz_p.h givaro.h givaro-zpz.h givaro-gfq.h lidia-gfq.h lidia.h

EXTRA_DIST=				\
	abstract.dxx			\
	fields.dxx			\
	parametric.dxx			\
	unparametric.dxx		\
	wrappers.dxx			\
	$(EXTRA_HEADERS)
