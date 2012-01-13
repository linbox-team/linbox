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

#INCLUDES=-I$(top_srcdir)
AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/linbox

#if LinBoxSrcOnly
#AM_CPPFLAGS=-g -DLinBoxSrcOnly $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
#LINBOX_LIBS=
#else
#AM_CPPFLAGS=-g $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
#LINBOX_LIBS= $(top_builddir)/linbox/liblinbox.la
#endif

# SAFER_FLAGS=--std=c++98 -g -Wall -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter
#  SAFER_FLAGS=-g -Wall -Wextra -Wno-unused-parameter

# AM_CXXFLAGS = $(SAFER_FLAGS) -DLinBoxSrcOnly -DDEBUG -O0
AM_CXXFLAGS = @DEFAULT_CFLAGS@ -DLinBoxSrcOnly
AM_CPPFLAGS += $(GMP_CFLAGS) $(NTL_CFLAGS) $(FFLAS_FFPACK_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(IML_CFLAGS)
#LDADD = $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS)  "-Wl,-enable-auto-import"
LDADD = $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS) $(MPFR_LIBS) $(FPLLL_LIBS) $(LDFLAGS) $(IML_LIBS)
#  AM_LDFLAGS=-static

#EX_LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS)
#LDADD = $(LINBOX_LIBS) $(EX_LDADD)

SUBDIRS=fields

EXAMPLES=rank det minpoly valence solve dot-product echelon sparseelimdet sparseelimrank checksolve doubledet smithvalence
# EXAMPLES+=nulp yabla
GIVARONTL_EXAMPLES=smith charpoly graph-charpoly
if LINBOX_HAVE_NTL
if LINBOX_HAVE_GIVARO
USE_GIVARONTL_EXAMPLES = $(GIVARONTL_EXAMPLES)
endif
endif


EXTRA_PROGRAMS= $(EXAMPLES) $(USE_GIVARONTL_EXAMPLES)

CLEANFILES=$(EXTRA_PROGRAMS)

EXTRA_DIST = examples.doxy

examples: ${EXTRA_PROGRAMS}

dot_product_SOURCES    = dot-product.C
doubledet_SOURCES      = doubledet.C
#linbox_program_SOURCES = \
#	linbox-program.C
checksolve_SOURCES     = checksolve.C
charpoly_SOURCES       = charpoly.C
graph_charpoly_SOURCES = graph-charpoly.C
det_SOURCES            = det.C
rank_SOURCES           = rank.C
smith_SOURCES          = smith.C
minpoly_SOURCES        = minpoly.C
valence_SOURCES        = valence.C
solve_SOURCES          = solve.C
echelon_SOURCES        = echelon.C
smithvalence_SOURCES   = smithvalence.C
sparseelimdet_SOURCES  = sparseelimdet.C
sparseelimrank_SOURCES = sparseelimrank.C
LINBOX=@prefix@

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=@bindir@


# for compilation of new examples
%:%.C
	$(CXX) $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDADD)

%:%.cpp
	$(CXX)  $(CXXFLAGS) $(AM_CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.cpp -o $@ $(LDADD)
