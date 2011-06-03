# Copyright (c) 2011 the LinBox group
# Brice Boyer <bboyer@imag.fr>
# This file is part of LinBox
# see COPYING for licence

AM_CPPFLAGS=-I$(top_srcdir) -I$(top_srcdir)/linbox

#  SAFER_FLAGS=-g -Wall -Wextra -Wno-unused-parameter

AM_CXXFLAGS= @DEFAULT_CFLAGS@ -DLinBoxSrcOnly
AM_CPPFLAGS+= $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS) $(FFLAFLAS_CFLAGS)
LDADD= $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS)
AM_LDFLAGS=-static $(LDFLAGS)

BENCH_BASIC=               \
		benchmark-fgemm    \
		benchmark-crafixed \
		benchmark-ftrXm
TODO=
		benchmark-matmul   \
		benchmark-spmv     \
		benchmark-fields

#  BENCH_ALGOS=               \
TODO=
		benchmark-solve    \
		benchmark-rank     \
		benchmark-det      \
		benchmark-nullspace

#  BENCH_FORMS=               \
TODO=
		benchmark-lu
		benchmark-echelon  \
		benchmark-hermite  \
		benchmark-smith

EXTRA_PROGRAMS= $(BENCH_BASIC)

EXTRA_DIST= benchmark.h benchmark.doxy

benchmarks: ${EXTRA_PROGRAMS}

### BASE BENCHMARK ###

benchmark_fgemm_SOURCES          = benchmark-fgemm.C
benchmark_ftrXm_SOURCES          = benchmark-ftrXm.C

benchmark_crafixed_SOURCES       = benchmark-crafixed.C
benchmark_crafixed_CPPFLAGS      = $(IML_CFLAGS) $(AM_CPPFLAGS)
benchmark_crafixed_LDADD         = $(IML_LIBS) $(LDADD)
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
