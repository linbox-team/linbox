# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence
pkgincludesubdir=$(pkgincludedir)/field

SUBDIRS=NTL Givaro FFPACK

BASIC_HDRS =            \
    field-interface.h   \
    field-traits.h      \
    rebind.h            \
    abstract.h          \
    archetype.h         \
    envelope.h          \
    unparametric.h      \
    param-fuzzy.h       \
    gmp-rational.h      \
    modular.h           \
    modular-balanced.h  \
    gf2.h               \
    gf2.inl             \
    hom.h               \
    map.h               \
    PIR-modular-int32.h \
    local2_32.h         \
    PID-integer.h       \
    PID-double.h        \
    block-ring.h        \
    multimod-field.h

NTL_HDRS =              \
    ntl.h               \
    PID-ntl-zz_p.h     \
    PIR-ntl-ZZ_p.h

GIVARO_HDRS =           \
    givaro.h

LIDIA_HDRS =            \
    lidia.h             \
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

pkgincludesub_HEADERS =     \
    $(BASIC_HDRS)           \
    $(USE_NTL_HDRS)         \
    $(USE_GIVARO_HDRS)      \
    $(USE_LIDIA_HDRS)

EXTRA_DIST =            \
    $(NTL_HDRS)         \
    $(GIVARO_HDRS)      \
    $(LIDIA_HDRS)       \
    field.doxy          \
    modular.doxy        \
    integers.doxy       \
    lidia.doxy          \
    ntl.doxy
