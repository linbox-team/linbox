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
# License along with this library; if not, write to the Free Software Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# ========LICENCE========



INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) $(ATLAS_CFLAGS) $(MAPLE_CFLAGS) $(GIVARO_CFLAGS)

EXTRA_DIST = lbmaple.C lbmaple.mpl.head lbmaple.mpl.tail Makefile.in.1 Makefile.in.2 demonstration.mws

lib_LTLIBRARIES=liblbmaple.la liblbmapleffpack.la

liblbmaple_la_SOURCES=lbmaple.C
liblbmaple_la_LDFLAGS=$(top_builddir)/linbox/liblinbox.la $(GMP_LIBS) $(NTL_LIBS) $(MAPLE_LIBS)

liblbmapleffpack_la_SOURCES=lbmaple-ffpack.C
liblbmapleffpack_la_LDFLAGS=$(top_builddir)/linbox/liblinbox.la $(ATLAS_LIBS) $(GMP_LIBS) $(NTL_LIBS) $(MAPLE_LIBS) $(GIVARO_LIBS)

install-data-local:
	$(INSTALL_DATA) lbmaple.mpl $(libdir)/lbmaple.mpl
