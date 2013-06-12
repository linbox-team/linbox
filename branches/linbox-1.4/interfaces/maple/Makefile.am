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



if LINBOX_HAVE_MAPLE
pkgincludesubdir=$(includedir)/interfaces/maple/new

AM_CPPFLAGS =-I$(top_srcdir) -I. -I../../linbox

AM_CXXFLAGS = @DEFAULT_CFLAGS@ -DDISABLE_COMMENTATOR -I../driver $(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(FFLAS_FFPACK_CFLAGS) $(MAPLE_CFLAGS) $(GIVARO_CFLAGS)

lib_LTLIBRARIES=liblbmaple.la

liblbmaple_la_SOURCES= lb-maple.C

liblbmaple_la_LDFLAGS=$(GIVARO_LIBS) $(GMP_LIBS) $(NTL_LIBS) $(BLAS_LIBS) $(MAPLE_LIBS) $(top_srcdir)/interfaces/driver/liblbdriver.la -Wl,-zmuldefs $(LDFLAGS)

pkginclude_HEADERS =lb-maple.h lb-maple-utilities.h


install-data-local:
	sh lb-maple-path.sh "$(libdir)/liblbmaple.so"
	$(INSTALL_DATA) lb-maple.mpl $(libdir)/lb-maple.mpl

install-data-hook:
	cd $(libdir) && \
	$(MAPLE_HOME)/bin/maple lb-maple.mpl 2>&1 > /dev/null
endif
