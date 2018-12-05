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

SUBDIRS=ntl modular

BASIC_HDRS =            \
    ring-interface.h    \
    polynomial-ring.h \
    modular.h           \
    local-pir-modular.h \
    pir-modular-int32.h \
    local2_32.h         \
    abstract.h          \
    archetype.h         \
    envelope.h          \
    polynomial-local-x.h

NTL_HDRS =              \
    ntl.h               \
    pid-ntl-zz_p.h     \
    pir-ntl-zz_p.h

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
    ntl.doxy
