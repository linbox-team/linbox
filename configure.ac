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

AC_INIT([LinBox], [1.5.2],[linbox-use@googlegroups.com],[linbox],
		[http://www.linalg.org/])

AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([1.10 -Wall -Wno-portability foreign])
AC_CONFIG_HEADERS([config.h])
AX_PREFIX_CONFIG_H(linbox/config.h, __LINBOX)

AM_PROG_AR

AM_ACLOCAL_INCLUDE(macros)

AC_LANG([C++])

echo "-----------------------------------------------"

AC_DEBUG
AC_PROFILE
AC_WARNINGS

echo "-----------------------------------------------"

# CFLAGS=${CFLAGS:-$DEFAULT_CFLAGS}
# CXXFLAGS=${CXXFLAGS:-$DEFAULT_CXXFLAGS}

######################################################
# Try and pass different flags according to compiler #
######################################################


# disable default -g -O2 CXXFLAGS
: ${CXXFLAGS=""}

#set CXX
AC_PROG_CXX

AC_COMPILER_NAME

STDFLAG=${CXXFLAGS}

CXXFLAGS=""

# We need a C++11 compiler now - AB 2014-12-12
# clang-3.8 does not support __float128 without explicitly passing it -std=c++11
if test "x${CCNAM}" = "xclang38" ; then
   AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
else
   AX_CXX_COMPILE_STDCXX_11([ext],[mandatory])
fi

AC_SUBST([DEFAULT_CFLAGS])
AC_SUBST([DEBUG_CFLAGS])
AC_SUBST([TESTS_CFLAGS])

TESTS_CFLAGS="-O2"
DEBUG_CFLAGS="-g"
DEFAULT_CFLAGS=""
WARN_CFLAGS="-Wall"

#TODO use -fast for icc, -ipa for eko...
if test "x$DBG" = "xyes" ; then
	DEFAULT_CFLAGS="-O0 ${DEFAULT_CFLAGS} " #those are CXXFLAGS
	DEBUG_CFLAGS="${DEBUG_CFLAGS} -DDEBUG -D_LB_DEBUG"
else
	DEFAULT_CFLAGS="-O2 ${DEFAULT_CFLAGS} "
	DEBUG_CFLAGS="${DEBUG_CFLAGS} -DNDEBUG -U_LB_DEBUG"
fi

if test "x$PROF" = "xyes" ; then
	DEFAULT_CFLAGS="${DEFAULT_CFLAGS} -pg"
fi

if test "x$WARN" = "xyes" -o "x$WARN" = "xfull" ; then
	case x${CCNAM} in
		xicc)
		WARN_CFLAGS="${WARN_CFLAGS} -Wcheck"
		WARN_CFLAGS="${WARN_CFLAGS} -Wall -Wno-unused-parameter -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wwrite-strings -Wno-long-long"
		WARN_CFLAGS="${WARN_CFLAGS} -Wextra -ansi"
			;;
		xeko)
			WARN_CFLAGS="${WARN_CFLAGS} -Wno-unused-parameter"
			;;
		xgcc|xgcc44)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -Wno-variadic-macros -Wno-vla"
			fi
			if test "x${HAVE_CXX11}" = "x0" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -ansi"
			fi
			;;
		xgcc48)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -Wno-variadic-macros -Wno-vla"
				# WARN_CFLAGS="${WARN_CFLAGS} -fsanitize=address"
			fi
			if test "x${HAVE_CXX11}" = "x0" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -ansi"
			fi

			;;

		xclang)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long -Wno-vla-extension -Wno-variadic-macros"
				WARN_CFLAGS="${WARN_CFLAGS} -D__STRICT_ANSI__"
			fi
			;;
		xclang31)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long  -g -Wno-vla-extension -Wno-variadic-macros"
				WARN_CFLAGS="${WARN_CFLAGS} -D__STRICT_ANSI__"
				# WARN_CFLAGS="${WARN_CFLAGS} -fsanitize=address"
			fi
			;;

		*)
			echo
			echo "*******************************************************"
			echo "unsupported compiler ($CCNAM). Please file a bug."
			echo "*******************************************************"
			echo
			WARN_CFLAGS="${WARN_CFLAGS}"
	esac
fi


DEFAULT_CFLAGS="${DEFAULT_CFLAGS} ${WARN_CFLAGS} ${DEBUG_CFLAGS}"
TESTS_CFLAGS="${TESTS_CFLAGS} ${WARN_CFLAGS} ${DEBUG_CFLAGS}"


AC_HEADER_STDC
AC_PROG_LIBTOOL
AC_PROG_SED

# newer libtool...
LT_PREREQ([2.4.2])
LT_INIT([win32-dll])


# needed for building interfaces as shared libs on Windows
AC_CANONICAL_HOST
case "${host_os}" in
    cygwin*|mingw*)
        LINBOX_LDFLAGS=-no-undefined
        ;;
    *)
        LINBOX_LDFLAGS=""
        ;;
esac
AC_SUBST([LINBOX_LDFLAGS])


echo "-----------------------------------------------"
echo "        START  LINBOX CONFIG                   "
echo "-----------------------------------------------"

# Checks which SIMD optimization flags to use and set SIMD_CFLAGS accordingly
INSTR_SET
AC_SUBST(SIMD_FLAGS)

dnl gcc-4.9.2 bug See https://trac.sagemath.org/ticket/17635#comment:178
AS_IF([ test  "x$CCNAM" = "xgcc492" ],[CXXFLAGS="${CXXFLAGS} -fpermissive"],[])


echo "-----------------------------------------------"
LB_OMP
echo "-----------------------------------------------"
AC_SUBST([PARFLAGS],['${OMPFLAGS}'])
AC_SUBST([PARLIBS],['${OMPFLAGS}'])

# Machine characteristics

AC_CHECK_SIZEOF(char, 8)
AC_CHECK_SIZEOF(short, 16)
AC_CHECK_SIZEOF(int, 32)
AC_CHECK_SIZEOF(long, 32)
AC_CHECK_SIZEOF(long long, 64)
AC_CHECK_SIZEOF(__int64, 64)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h sys/time.h stdint.h pthread.h])


# check endianness of the architecture
AC_C_BIGENDIAN(
  [AC_DEFINE(HAVE_BIG_ENDIAN, 1, [Define that architecture uses big endian storage])],
  [AC_DEFINE(HAVE_LITTLE_ENDIAN, 1, [Define that architecture uses little endian storage])],
  [])

echo "-----------------------------------------------"
# Feature checks
LB_MISC
LB_DRIVER


echo "-----------------------------------------------"

# Now getting GMP and Givaro from FFLAS-FFPACK - AB 2014-12-10
PKG_CHECK_MODULES([FFLAS_FFPACK], [fflas-ffpack >= 2.3.0])

dnl LB_CHECK_FFLAS_FFPACK(,,[
dnl echo ''
dnl echo '*******************************************************************************'
dnl echo ' ERROR: Fflas-Ffpack not found!'
dnl echo
dnl echo ' Fflas-Ffpack routines are required for this library to compile. Please'
dnl echo ' make sure they are installed and specify its location with the option'
dnl echo ' --with-fflas-ffpack=<lib> when running configure.'
dnl echo ' Also make sure your compiler supports cxx-11...'
dnl echo '*******************************************************************************'
dnl exit 1
dnl ])

LB_CHECK_LAPACK

# AC_MSG_RESULT($BLAS_LIBS)
# AC_MSG_RESULT($BLAS_CFLAGS)

LB_CHECK_IML
LB_CHECK_M4RI
LB_CHECK_M4RIE
LB_CHECK_MPFR
# we nedd mpfr next :
LB_CHECK_FPLLL
LB_CHECK_FLINT

LB_CHECK_NTL

LB_CHECK_SACLIB
LB_CHECK_MAPLE
LB_CHECK_SAGE

LB_CHECK_EXPAT

LB_BENCH
LB_CHECK_XML
# LB_OPT
LB_DOC

LB_CHECK_OCL


if test ! -d ./benchmarks/data ; then
	echo "Creating data dir in benchmark" ;
	mkdir ./benchmarks/data ;
fi

DEPS_CFLAGS="${FFLAS_FFPACK_CFLAGS} ${NTL_CFLAGS} ${MPFR_CFLAGS} ${FPLLL_CFLAGS} ${IML_CFLAGS} ${FLINT_CFLAGS}"
DEPS_LIBS=" ${NTL_LIBS} ${MPFR_LIBS} ${FPLLL_LIBS} ${IML_LIBS} ${FLINT_LIBS} ${OCL_LIBS} ${FFLAS_FFPACK_LIBS} ${XML_LIBS}"

CXXFLAGS="${CXXFLAGS} ${STDFLAG}"
CXXFLAGS="${CXXFLAGS} ${SIMD_CFLAGS}"

AC_SUBST(CXXFLAGS)
AC_SUBST(DEPS_CFLAGS)
AC_SUBST(DEPS_LIBS)

echo "-----------------------------------------------"

dnl MYLIBS="${FFLAS_FFPACK_LIBS} ${NTL_LIBS} ${MPFR_LIBS} ${FPLLL_LIBS} ${IML_LIBS} ${FLINT_LIBS} ${LDFLAGS} ${OCL_LIBS}"
dnl MYLIBS="${LIBS}"
dnl MYFLAGS="${FFLAS_FFPACK_CFLAGS} ${NTL_CFLAGS} ${MPFR_CFLAGS} ${FPLLL_CFLAGS} ${IML_CFLAGS} ${FLINT_CFLAGS} ${OPT_FLAGS}"

echo "CXXFLAGS       = ${CXXFLAGS}"
echo "LIBS           = ${LIBS}"
echo "DEPS_CFLAGS    = ${DEPS_CFLAGS}"
echo "DEPS_LIBS      = ${DEPS_LIBS}"
echo "DEFAULT CFLAGS = ${DEFAULT_CFLAGS}"
echo "TESTS CFLAGS   = ${TESTS_CFLAGS}"

echo "-----------------------------------------------"
echo "         END   LINBOX CONFIG                   "
echo "-----------------------------------------------"

AC_CONFIG_FILES([
Makefile
linbox-config
examples/Makefile
examples/data/Makefile
doc/Makefile
linbox/Makefile
linbox/algorithms/Makefile
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
interfaces/sage/Makefile
macros/Makefile
macros/CodeChunk/Makefile
benchmarks/Makefile
benchmarks/data/Makefile
benchmarks/matrix/Makefile
linbox.pc
])
# linbox/algorithms/IML/Makefile

# examples/Readme-make
# examples/field
# examples/solver
# examples/blackbox
# examples/data

AC_OUTPUT

echo "-----------------------------------------------"
