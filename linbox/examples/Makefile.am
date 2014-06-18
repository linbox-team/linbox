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

DEFAULT_INCLUDES=-I@includedir@

#  -DLinBoxSrcOnly
DEFCPPFLAGS  = @DEFAULT_CFLAGS@   -DDISABLE_COMMENTATOR
AM_CPPFLAGS = $(DEFCPPFLAGS) $(GMP_CFLAGS) $(NTL_CFLAGS) $(BLAS_CFLAGS) $(FFLAS_FFPACK_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(IML_CFLAGS) $(OPTFLAGS)
LDADD = $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS) $(MPFR_LIBS) $(FPLLL_LIBS) $(LDFLAGS) $(IML_LIBS)
#  AM_LDFLAGS=-static
AM_LDFLAGS=-L@libdir@ -llinbox
if LINBOX_HAVE_OCL
AM_CPPFLAGS += $(OCL_CFLAGS)
LDADD += $(OCL_LIBS)
endif



#EX_LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS)
#LDADD = $(LINBOX_LIBS) $(EX_LDADD)

#  SUBDIRS=fields solver data blackbox

EXAMPLES=rank det minpoly valence solve dot-product echelon sparseelimdet sparseelimrank checksolve doubledet smithvalence charpoly polysmith
# EXAMPLES+=nulp yabla
GIVARONTL_EXAMPLES=smith graph-charpoly
if LINBOX_HAVE_NTL
if LINBOX_HAVE_GIVARO
USE_GIVARONTL_EXAMPLES = $(GIVARONTL_EXAMPLES)
endif
endif


EXTRA_PROGRAMS= $(EXAMPLES) $(USE_GIVARONTL_EXAMPLES)

CLEANFILES=$(EXTRA_PROGRAMS)

EXTRA_DIST = examples.doxy smithvalence.h
#map-sparse.h map-sparse.inl

examples:  ${EXTRA_PROGRAMS}

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
polysmith_SOURCES      = poly-smith.C
LINBOX=@prefix@

LINBOX_BIN=@bindir@


# for compilation of new examples
%:%.C
	$(CXX) -I@includedir@ $(CXXFLAGS) $(AM_CPPFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ -L@libdir@ -llinbox $(AM_LDFLAGS) $(LDADD)

%:%.cpp
	$(CXX)  $(CXXFLAGS) $(AM_CPPFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.cpp -o $@ $(AM_LDFLAGS) $(LDADD)
