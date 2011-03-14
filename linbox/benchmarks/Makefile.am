# Copyright (c) 2011 the LinBox group
# This file is part of LinBox
# see COPYING for licence

#INCLUDES=-I$(top_srcdir)
INCLUDES=-I$(top_srcdir) -I$(top_srcdir)/linbox

# SAFER_FLAGS=--std=c++98 -g -Wall -Wextra -Wno-ignored-qualifiers -Wno-unused-parameter
SAFER_FLAGS=-g -Wall -Wextra -Wno-unused-parameter

# AM_CXXFLAGS= $(SAFER_FLAGS) -DLinBoxSrcOnly -DDEBUG -O0
AM_CXXFLAGS= $(SAFER_FLAGS) -DLinBoxSrcOnly -O2
AM_CPPFLAGS= $(GMP_CFLAGS) $(NTL_CFLAGS) $(GIVARO_CFLAGS) $(LIDIA_CFLAGS)
#LDADD= $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS)  "-Wl,-enable-auto-import"
LDADD= $(NTL_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS) $(GMP_LIBS)

#EX_LDADD = $(NTL_LIBS) $(GMP_LIBS) $(GIVARO_LIBS) $(LIDIA_LIBS) $(BLAS_LIBS)
#LDADD = $(LINBOX_LIBS) $(EX_LDADD)

#  BENCHS=                    \
		#  benchmark-matmul   \
		#  benchmark-fgemm

EXTRA_PROGRAMS= $(BENCHS)

#  EXTRA_DIST= benchmark.h benchmark.doxy

benchmarks: ${EXTRA_PROGRAMS}

#  benchmark_matmul_SOURCES         = benchmark-matmul.C
#  benchmark_fgemm_SOURCES          = benchmark-fgemm.C

LINBOX=@prefix@

# no /bin if LinBox is uncompiled
#LINBOX_BIN=$(LINBOX)/bin
LINBOX_BIN=@bindir@


# for compilation of new examples
.C:
	$(CXX) $(SAFER_FLAGS) -O2 $(CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.C -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(top_srcdir)/linbox/.libs/liblinbox.a

.cpp:
	$(CXX) $(SAFER_FLAGS) $(CXXFLAGS) $(OPTFLAGS) ${INCLUDES} $(AM_CPPFLAGS) $*.cpp -o $@ $(LDFLAGS) $(LDADD) $(LOADLIBES) $(top_srcdir)/linbox/.libs/liblinbox.a
