# Copyright (c) 2010 the LinBox group
# ========LICENCE========
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
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========



pkgincludesubdir=$(pkgincludedir)/ring

SUBDIRS=NTL Modular

BASIC_HDRS =            \
    ring-interface.h    \
    givaro-poly.h       \
    givaro-polynomial.h \
    modular.h           \
    PIR-modular-int32.h \
    local2_32.h         \
    PID-integer.h       \
    abstract.h          \
    archetype.h         \
    envelope.h          \
    poweroftwomodular.h

NTL_HDRS =              \
    ntl.h               \
    PID-ntl-zz_p.h     \
    PIR-ntl-ZZ_p.h

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

pkgincludesub_HEADERS =     \
    $(BASIC_HDRS)           \
    $(USE_NTL_HDRS)

EXTRA_DIST =            \
    $(NTL_HDRS)         \
    ring.doxy           \
    polynomial.doxy		\
    integers.doxy       \
    ntl.doxy
