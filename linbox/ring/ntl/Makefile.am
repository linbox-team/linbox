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



pkgincludesubdir=$(pkgincludedir)/ring/ntl


NTL_HDRS =              \
    ntl-rr.h            \
    ntl-zz_p.h          \
    ntl-zz_pe.h         \
    ntl-zz_px.h         \
    ntl-zz.h            \
    ntl-lzz_p.h         \
    ntl-lzz_pe.h        \
    ntl-lzz_px.h        \
    ntl-lzz_pex.h       \
    ntl-gf2e.h


if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif


pkgincludesub_HEADERS =     \
    $(USE_NTL_HDRS)

#  EXTRA_DIST =            \
    #  ntl.doxy
