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



pkgincludesubdir=$(pkgincludedir)/field/NTL


NTL_HDRS =              \
    ntl-RR.h            \
    ntl-ZZ_p.h          \
    ntl-ZZ_pE.h         \
    ntl-ZZ_pX.h         \
    ntl-ZZ.h            \
    ntl-lzz_p.h         \
    ntl-lzz_pE.h        \
    ntl-lzz_pX.h        \
    ntl-lzz_pEX.h       \
    ntl-GF2E.h


if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif


pkgincludesub_HEADERS =     \
    $(USE_NTL_HDRS)

#  EXTRA_DIST =            \
    #  ntl.doxy
