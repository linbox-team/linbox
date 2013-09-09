# Copyright (c) 2011 the LinBox group
# Brice Boyer <bboyer@imag.fr>
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

SUBDIRS=data matrix

pkgincludesubdir=$(pkgincludedir)/benchmarks

#  noinst_LTLIBRARIES=libbenchmarks.la
#  libbenchmarks_la_SOURCES= benchmark.C
#  libbenchmarks_la_LDFLAGS=$(LDFLAGS)

#  libbenchmarks_la_LIBADD=$(top_srcdir)/linbox/.libs/liblinbox.la


AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/linbox
AM_CPPFLAGS += $(GMP_CFLAGS) $(NTL_CFLAGS)  $(FFLAS_FFPACK_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(IML_CFLAGS)
AM_CXXFLAGS = @DEFAULT_CFLAGS@ -DDISABLE_COMMENTATOR -DLinBoxSrcOnly
LDADD = $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS) $(MPFR_LIBS) $(FPLLL_LIBS) $(LDFLAGS) $(IML_LIBS)
#  AM_LDFLAGS=-static $(LDFLAGS)
#  AM_LDFLAGS = .libs/libbenchmarks.la
#  AM_LDFLAGS = .libs/benchmark.o



BENCH_BASIC=               \
		benchmark-example

FAILS=    \
		benchmark-ftrXm \
		benchmark-ftrXm \
		benchmark-crafixed

TODO= \
		benchmark-matmul   \
		benchmark-spmv     \
		benchmark-fields

#  BENCH_ALGOS=               \
TODO= \
		benchmark-solve    \
		benchmark-rank     \
		benchmark-det      \
		benchmark-nullspace

#  BENCH_FORMS=               \
TODO= \
		benchmark-lu
		benchmark-echelon  \
		benchmark-hermite  \
		benchmark-smith

EXTRA_PROGRAMS= $(BENCH_BASIC)

pkginclude_HEADERS = \
		     benchmark.h \
		     benchmark.C \
		     benchmark.inl

EXTRA_DIST= \
	    benchmark.doxy

benchmarks: ${EXTRA_PROGRAMS}

### BASE BENCHMARK ###

#  benchmark_fgemm_SOURCES         = benchmark-fgemm.C

#  benchmark_ftrXm_SOURCES         = benchmark-ftrXm.C

#  benchmark_crafixed_SOURCES      = benchmark-crafixed.C

benchmark_example_SOURCES       = benchmark-example.C

#  benchmark_matmul_SOURCES         = benchmark-matmul.C
#  benchmark_spmv_SOURCES           = benchmark-spmv.C
#  benchmark_fields_SOURCES         = benchmark-fields.C

### BENCHMARK ALGOS and SOLUTIONS ###
#  benchmark_solve_SOURCES          = benchmark-solve.C
#  benchmark_rank_SOURCES           = benchmark-rank.C
#  benchmark_det_SOURCES            = benchmark-det.C
#  benchmark_nullspace_SOURCES      = benchmark-nullspace.C


### BENCHMARK MATRIX FACTORISATIONS ###
#  benchmark_lu_SOURCES             = benchmark-lu.C
#  benchmark_echelon_SOURCES        = benchmark-echelon.C
#  benchmark_hermite_SOURCES        = benchmark-hermite.C
#  benchmark_smith_SOURCES          = benchmark-smith.C


LINBOX=@prefix@

LINBOX_BIN=@bindir@

# for compilation of new benchmarks
%:%.C
	$(CXX) $(AM_CXXFLAGS) $(CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(top_srcdir)/linbox/.libs/liblinbox.a

%:%.cpp
	$(CXX)  $(AM_CXXFLAGS) $(CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.cpp -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(top_srcdir)/linbox/.libs/liblinbox.a
