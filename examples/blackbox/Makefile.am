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

EXTRA_PROGRAMS = smith example ex-mat0 load-det load-minpoly load-symmetrize random-rank
#zeroone-rank

AM_CPPFLAGS=-DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) -DEXAMPLE_DATADIR=\"$(srcdir)/data\"
LDADD = $(GMP_LIBS) $(NTL_LIBS) -lgivaro -lclapack -lcblas -latlas  $(LDFLAGS)

example_SOURCES  = example.C
smith_SOURCES  = smith.C
load_det_SOURCES = load-det.C
load_minpoly_SOURCES = load-minpoly.C
load_symmetrize_SOURCES = load-symmetrize.C
random_rank_SOURCES = random-rank.C
#zeroone_rank_SOURCES = zeroone-rank.C
ex_mat0_SOURCES  = ex-mat0.C

# general

%:%.o
	$(CXXLINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)


