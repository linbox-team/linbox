includedir=$(prefix)/include/LinBox

lib_LIBRARIES=liblinbox_givaro.a

include_HEADERS=\
	givcompgmp.h		giverror.h		\
	givinteger.h		givintfactor.h		\
	givintnumtheo.h		givintprime.h		\
	givintrns.h		givintrsa.h		\
	givinteger.inl		givintfactor.inl	\
	givintnumtheo.inl	givintrns_convert.inl	\
	givintrns_cstor.inl	givintrsa.inl

liblinbox_givaro_a_SOURCES=\
	giverror.C	givintadd.C	\
	givintcompare.C	givintcstor.C	\
	givintdiv.C	givinteger.C	\
	givintgcd.C	givintio.C	\
	givintmisc.C	givintmod.C	\
	givintmul.C	givintnumtheo.C	\
	givintpow.C	givintprime.C	\
	givintrsa.C	givintsub.C

INCLUDES=\
	$(GMP_CFLAGS)