# Copyright (c) 2010 the LinBox group
# This file is part of LinBox (and probably stolen here and there)
#  ========LICENCE========
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

AC_PREREQ([2.68])

AC_INIT([LinBox], [1.7.0],[linbox-use@googlegroups.com],[linbox],
		[http://www.linalg.org/])

AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign subdir-objects])
AX_PREFIX_CONFIG_H(linbox/config.h, __LINBOX)

# Set CXXFLAGS to an empty string if not defined, and save it. We must do it
# before calling AC_PROG_CXX that sets it to '-O2 -g' if not defined
: ${CXXFLAGS=""}
ORIGINAL_CXXFLAGS="${CXXFLAGS}"

# We set the language to C++
AC_LANG([C++])
AC_PROG_CXX

#
AM_PROG_AR

# Libtool
AC_PROG_LIBTOOL
LT_PREREQ([2.4.3])
LT_INIT([win32-dll])

# Look for headers
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])

#
AC_PROG_EGREP
AC_PROG_SED

#################################################
AS_BOX([          START LINBOX CONFIG          ])
#################################################

AC_COMPILER_NAME

# We need a C++11 compiler now - AB 2014-12-12
# clang-3.8 does not support __float128 without explicitly passing it -std=c++11
AS_IF([test "x${CCNAM}" = "xclang38"],
      [AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])],
      [AX_CXX_COMPILE_STDCXX_11([ext],[mandatory])]
      )
REQUIRED_FLAGS="${CXX11FLAGS}"
# Add the c++11 flags for the configure compilations as clang needs it to work properly with a glibc++ compiled with gcc
CXXFLAGS="${CXX11FLAGS} ${CXXFLAGS}"

AS_ECHO([---------------------------------------])

# Set OPTIM_FLAGS, DEBUG_FLAGS depending on compiler and command line arguments
SET_FLAGS

# Append -march=native or -mcpu=native (if recognized by the compiler) to
# OPTIM_FLAGS if not present in CXXFLAGS and not cross-compiling and
# --without-archnative is not set
ARCH_FLAGS

# Append -mfpmath=sse to OPTIM_FLAGS on i386 and i686 architecture with SSE
FPMATH_FLAGS

AS_ECHO([---------------------------------------])
# Machine characteristics

# Size of some types
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(__int64_t)

# Looking for int128
AC_CHECK_TYPE([__int128_t], [AC_DEFINE(HAVE_INT128, 1, [Define that compiler allows int128_t types])])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

AS_ECHO([---------------------------------------])

# Feature checks
LB_MISC
LB_DRIVER

# Looking for FFLAS-FFPACK.
# We get the flags for OpenMP, Givaro and GMP at the same time
PKG_CHECK_MODULES([FFLAS_FFPACK], [fflas-ffpack >= 2.5.0])
AC_SUBST(FFLAS_FFPACK_CFLAGS)
AC_SUBST(FFLAS_FFPACK_LIBS)

AS_ECHO([---------------------------------------])

LB_CHECK_LAPACK

LB_CHECK_NTL
LB_CHECK_MPFR
LB_CHECK_FPLLL
LB_CHECK_IML
LB_CHECK_FLINT

LB_CHECK_SACLIB
LB_CHECK_MAPLE

AS_ECHO([---------------------------------------])

LB_BENCH
LB_CHECK_XML
# LB_OPT
LB_DOC

AS_ECHO([---------------------------------------])

LB_CHECK_OCL
LB_CHECK_MPI

AS_ECHO([---------------------------------------])

# needed for building interfaces as shared libs on Windows
AC_CANONICAL_HOST
AS_CASE([$host_os], [cygwin*|mingw*], [ LINBOX_LDFLAGS="-no-undefined" ],
                                      [ LINBOX_LDFLAGS="" ])
AC_SUBST([LINBOX_LDFLAGS])

AS_ECHO([---------------------------------------])

CXXFLAGS="${ORIGINAL_CXXFLAGS}"
LINBOX_CXXFLAGS="${REQUIRED_FLAGS} ${OPTIM_FLAGS} ${DEBUG_FLAGS}"
AS_ECHO(["LINBOX_CXXFLAGS = ${LINBOX_CXXFLAGS}"])

AC_SUBST(LINBOX_CXXFLAGS)
AC_SUBST(REQUIRED_FLAGS)

LINBOX_DEPS_CFLAGS="${NTL_CFLAGS} ${MPFR_CFLAGS} ${FPLLL_CFLAGS} ${IML_CFLAGS} ${FLINT_CFLAGS} ${OCL_CFLAGS}"
LINBOX_DEPS_LIBS="${NTL_LIBS} ${MPFR_LIBS} ${FPLLL_LIBS} ${IML_LIBS} ${FLINT_LIBS} ${OCL_LIBS}"

AC_SUBST(LINBOX_DEPS_CFLAGS)
AC_SUBST(LINBOX_DEPS_LIBS)

#################################################
AS_BOX([           END LINBOX CONFIG           ])
#################################################

# symlinks to matrix files for the testsuite
AC_CONFIG_LINKS([tests/data/sms.matrix:tests/data/sms.matrix])
AC_CONFIG_LINKS([tests/data/30_30_27.sms:tests/data/30_30_27.sms])
AC_CONFIG_LINKS([tests/data/fib25.sms:tests/data/fib25.sms])
AC_CONFIG_LINKS([tests/data/rk9_7_10.sms:tests/data/rk9_7_10.sms])

AC_CONFIG_FILES([
Makefile
linbox-config
examples/Makefile
examples/data/Makefile
doc/Makefile
linbox/Makefile
linbox/algorithms/Makefile
linbox/algorithms/dixon-solver/Makefile
linbox/algorithms/gauss/Makefile
linbox/algorithms/matrix-blas3/Makefile
linbox/algorithms/opencl-kernels/Makefile
linbox/algorithms/polynomial-matrix/Makefile
linbox/blackbox/Makefile
linbox/element/Makefile
linbox/field/Makefile
linbox/matrix/Makefile
linbox/matrix/sparsematrix/Makefile
linbox/matrix/densematrix/Makefile
linbox/matrix/matrixdomain/Makefile
linbox/matrix/sliced3/Makefile
linbox/randiter/Makefile
linbox/ring/Makefile
linbox/ring/ntl/Makefile
linbox/ring/modular/Makefile
linbox/solutions/Makefile
linbox/solutions/solve/Makefile
linbox/solutions/echelon/Makefile
linbox/util/Makefile
linbox/util/formats/Makefile
linbox/vector/Makefile
linbox/polynomial/Makefile
tests/Makefile
tests/data/Makefile
tests/matrix/Makefile
interfaces/Makefile
interfaces/driver/Makefile
interfaces/maple/Makefile
interfaces/kaapi/Makefile
benchmarks/Makefile
benchmarks/data/Makefile
benchmarks/matrix/Makefile
linbox.pc
])

# examples/Readme-make
# examples/field
# examples/solver
# examples/blackbox
# examples/data

AC_OUTPUT
