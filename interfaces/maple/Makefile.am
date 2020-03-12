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


AM_CPPFLAGS=-I$(top_srcdir) -I. -I../../linbox -I../driver -DDISABLE_COMMENTATOR
AM_CXXFLAGS=$(LINBOX_CXXFLAGS) $(FFLAS_FFPACK_CFLAGS) $(LINBOX_DEPS_CFLAGS) $(MAPLE_CFLAGS)

lib_LTLIBRARIES=liblbmaple.la

liblbmaple_la_SOURCES= lb-maple.C

liblbmaple_la_LDFLAGS=$(LINBOX_DEPS_LIBS) $(FFLAS_FFPACK_LIBS) $(LINBOX_LDFLAGS) $(MAPLE_LIBS) $(top_srcdir)/interfaces/driver/liblbdriver.la -Wl,-zmuldefs

pkginclude_HEADERS =lb-maple.h lb-maple-utilities.h


install-data-local:
	sh lb-maple-path.sh "$(libdir)/liblbmaple.so"
	$(INSTALL_DATA) lb-maple.mpl $(libdir)/lb-maple.mpl

install-data-hook:
	cd $(libdir) && \
	$(MAPLE_HOME)/bin/maple lb-maple.mpl 2>&1 > /dev/null
endif
