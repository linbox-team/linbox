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

AM_CPPFLAGS=-I$(top_srcdir) -DDISABLE_COMMENTATOR
AM_CXXFLAGS=$(LINBOX_CXXFLAGS) $(FFLAS_FFPACK_CFLAGS) $(LINBOX_DEPS_CFLAGS)
LDADD=$(top_builddir)/linbox/liblinbox.la $(LINBOX_DEPS_LIBS) $(FFLAS_FFPACK_LIBS) $(LINBOX_LDFLAGS)


#  SUBDIRS=fields solver data blackbox
SUBDIRS=data

EXAMPLES=rank det minpoly valence solve dot-product echelon sparseelimdet \
sparseelimrank checksolve doubledet smithvalence charpoly blassolve \
dixonsolve sparsesolverat densesolverat \
ratdet poweroftwo_ranks power_rank genprime smithsparse nullspacebasis_rat nullspacebasis matrices
#polysmith bench-fft bench-matpoly-mult
# EXAMPLES+=nulp yabla 
GIVARONTL_EXAMPLES=smith graph-charpoly
if LINBOX_HAVE_NTL
USE_GIVARONTL_EXAMPLES = $(GIVARONTL_EXAMPLES)
endif


EXTRA_PROGRAMS= $(EXAMPLES) $(USE_GIVARONTL_EXAMPLES)

CLEANFILES=$(EXTRA_PROGRAMS)

EXTRA_DIST = examples.doxy matrices.h matrices.C
#map-sparse.h map-sparse.inl

examples:  ${EXTRA_PROGRAMS}

dot_product_SOURCES    = dot-product.C
doubledet_SOURCES      = doubledet.C
checksolve_SOURCES     = checksolve.C
charpoly_SOURCES       = charpoly.C
graph_charpoly_SOURCES = graph-charpoly.C
det_SOURCES            = det.C
genprime_SOURCES       = genprime.C
rank_SOURCES           = rank.C
smith_SOURCES          = smith.C
smithsparse_SOURCES    = smithsparse.C
matrices_SOURCES       = matrices.C
minpoly_SOURCES        = minpoly.C
valence_SOURCES        = valence.C
solve_SOURCES          = solve.C
echelon_SOURCES        = echelon.C
smithvalence_SOURCES   = smithvalence.C
sparseelimdet_SOURCES  = sparseelimdet.C
sparseelimrank_SOURCES = sparseelimrank.C
dixonsolve_SOURCES  = dixonsolve.C
ratdet_SOURCES       = ratdet.C
sparsesolverat_SOURCES = sparsesolverat.C
densesolverat_SOURCES = densesolverat.C
blassolve_SOURCES      = blassolve.C
power_rank_SOURCES     = power_rank.C
poweroftwo_ranks_SOURCES=poweroftwo_ranks.C
nullspacebasis_rat_SOURCES=nullspacebasis_rat.C
nullspacebasis_SOURCES=nullspacebasis.C
#smithformlocal_SOURCES = smith-form-local.C
#polysmith_SOURCES      = poly-smith.C
#bench_matpoly_mult_SOURCES = bench-matpoly-mult.C

LINBOX=@prefix@

LINBOX_BIN=@bindir@


# for compilation of new examples
define comp_new_examp
	$(AM_V_CXX)$(CXXCOMPILE) -c -o $@.$(OBJEXT) $<
	$(AM_V_CXXLD)$(CXXLINK) $@.$(OBJEXT) $(LDADD)
endef	

%:%.C
	$(comp_new_examp)

%:%.cpp
	$(comp_new_examp)
