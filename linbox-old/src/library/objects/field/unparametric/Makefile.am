includedir=$(prefix)/include/LinBox

EXTRA_DIST=unparametric.dxx

include_HEADERS=\
	fuzzy.h				\
	gmp-rational-field.h		\
	gmp-rational-field.C		\
	modular.h

lib_LIBRARIES=liblinbox_gmp_rational_field.a

liblinbox_gmp_rational_field_a_SOURCES=gmp-rational-field.C

INCLUDES=\
	$(GMP_FLAGS)

