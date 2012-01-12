# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# ========LICENCE========
# This file is part of the library FFLAS-FFPACK.
#
# FFLAS-FFPACK is free software: you can redistribute it and/or modify
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

AM_CPPFLAGS=-I$(top_srcdir)

EXTRA_PROGRAMS =		\
	ex-fields 		\
	ex-fields-archetype	\
	modular-int

#AM_CPPFLAGS=$(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
#LDADD = $(top_builddir)/linbox/liblinbox.a $(GMP_LIBS) $(NTL_LIBS)
AM_CXXFLAGS = @DEFAULT_CFLAGS@ -DLinBoxSrcOnly
AM_CPPFLAGS+= -DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) \
			 -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(GMP_LIBS) $(NTL_LIBS) -lgivaro -lclapack -lcblas -latlas $(LDFLAGS)

ex_fields_SOURCES           = ex-fields.C
ex_fields_archetype_SOURCES = ex-fields-archetype.C
modular_int_SOURCES = modular-int.C
