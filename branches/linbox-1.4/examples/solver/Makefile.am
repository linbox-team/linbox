# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
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
#/

INCLUDES=-I$(top_srcdir)

AM_CPPFLAGS=-g $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(ATLAS_CFLAGS)
LDADD = $(top_builddir)/linbox/liblinbox.la $(NTL_LIBS) $(GIVARO_LIBS)  $(GMP_LIBS) $(LIDIA_LIBS) $(ATLAS_LIBS) $(LDFLAGS)

EXTRA_PROGRAMS=t-rdisolve

t_rdisolve_SOURCES = \
	t-rdisolve.C

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=$(bindir)

# for compilation of new examples
%:%.C
	$(CXX) $(CPPFLAGS) $*.C -o $*  $(shell $(LINBOX_BIN)/linbox-config --cflags --libs) $(LDADD)
