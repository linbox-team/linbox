# Copyright (c) 2010 the LinBox group
#  ========LICENCE========
# This file is part of the library LinBox.
#
# LinBox is free software: you can redistribute it and/or modify
# it under the terms of the  GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========



pkgincludesubdir=$(pkgincludedir)/field

SUBDIRS=NTL Givaro Modular

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
