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



pkgincludesubdir=$(pkgincludedir)/field/Givaro

GIVARO_HDRS =           \
    givaro-zpz.h        \
    givaro-montg.h      \
    givaro-zpz.inl      \
    givaro-gfq.h        \
    givaro-extension.h  \
    givaro-field.h      \
    givaro-field-gf2.h  \
    givaro-rational.h

if LINBOX_HAVE_GIVARO
USE_GIVARO_HDRS = $(GIVARO_HDRS)
endif

pkgincludesub_HEADERS =     \
    $(USE_GIVARO_HDRS)

#  EXTRA_DIST =            \
    #  givaro.doxy
