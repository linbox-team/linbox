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


AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/linbox   
LIBS=
AM_CXXFLAGS=@DEFAULT_CFLAGS@
pkgincludesubdir=$(pkgincludedir)/randiter
 
#AM_CPPFLAGS=-I$(top_srcdir)/linbox
#AM_CXXFLAGS = @DEFAULT_CFLAGS@ -DDISABLE_COMMENTATOR $(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS)   

BASIC_HDRS =            \
    abstract.h          \
    archetype.h         \
    envelope.h          \
    gmp-rational.h      \
    param-fuzzy.h       \
    modular-balanced.h  \
    modular-crooked.h   \
    nonzero.h           \
    mersenne-twister.h  \
    ntl-ZZ.h            \
    random-prime.h      \
    random-integer.h    \
    gmp-random-prime.h  \
    random-fftprime.h   \
    multimod-randomprime.h \
    givaro-poly.h

NTL_HDRS =

LIDIA_HDRS =            \
    lidia.h             \
    lidia-gfq.h

if LINBOX_HAVE_NTL
USE_NTL_HDRS = $(NTL_HDRS)
endif

if LINBOX_HAVE_LIDIA
USE_LIDIA_HDRS = $(LIDIA_HDRS)
endif

pkgincludesub_HEADERS =     \
    $(BASIC_HDRS)           \
    $(USE_NTL_HDRS)         \
    $(USE_LIDIA_HDRS)

EXTRA_DIST =                \
    $(NTL_HDRS)             \
    $(LIDIA_HDRS)           \
    randiter.doxy           \
    archetype.doxy          \
    randiter-wrappers.doxy

noinst_LTLIBRARIES = libranditer.la

libranditer_la_SOURCES =    \
    mersenne-twister.C

