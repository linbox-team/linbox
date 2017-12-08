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

# DEFAULT_INCLUDES=-I@includedir@

#  -DLinBoxSrcOnly
DEFCPPFLAGS  = @DEFAULT_CFLAGS@   -DDISABLE_COMMENTATOR
AM_CPPFLAGS = -I$(top_srcdir) $(DEFCPPFLAGS) $(DEPS_CFLAGS)
AM_CPPFLAGS += $(OPTFLAGS) 
LDADD = $(DEPS_LIBS) $(LDFLAGS)
LDADD += $(top_builddir)/linbox/liblinbox.la
LDADD += $(OPTLIBES)
#  AM_LDFLAGS=-static

#AM_LDFLAGS=-L@libdir@ -llinbox

# if LINBOX_HAVE_OCL
# AM_CPPFLAGS += $(OCL_CFLAGS)
# LDADD += $(OCL_LIBS)
# endif



#EX_LDADD = $(NTL_LIBS)
#LDADD = $(LINBOX_LIBS) $(EX_LDADD)

#  SUBDIRS=fields solver data blackbox
SUBDIRS=data

EXAMPLES=rank det minpoly valence solve dot-product echelon sparseelimdet sparseelimrank checksolve doubledet smithvalence charpoly blassolve solverat sparsesolverat poweroftwo_ranks power_rank#polysmith bench-fft bench-matpoly-mult
# EXAMPLES+=nulp yabla 
GIVARONTL_EXAMPLES=smith graph-charpoly
if LINBOX_HAVE_NTL
USE_GIVARONTL_EXAMPLES = $(GIVARONTL_EXAMPLES)
endif


EXTRA_PROGRAMS= $(EXAMPLES) $(USE_GIVARONTL_EXAMPLES)

CLEANFILES=$(EXTRA_PROGRAMS)

EXTRA_DIST = examples.doxy smithvalence.h test.sh
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
solverat_SOURCES       = solverat.C
sparsesolverat_SOURCES = sparsesolverat.C
blassolve_SOURCES      = blassolve.C
power_rank_SOURCES     = power_rank.C
poweroftwo_ranks_SOURCES=poweroftwo_ranks.C
#smithformlocal_SOURCES = smith-form-local.C
#polysmith_SOURCES      = poly-smith.C
#bench_fft_SOURCES       = bench-fft.C
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
