# Copyright (c) 2010 the LinBox group
# This file is part of LinBox (and probably stolen here and there)
#  ========LICENCE========
# This file is part of the library LinBox.
#
 LinBox is free software: you can redistribute it and/or modify
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

AC_PREREQ([2.64])

AC_INIT([LinBox], [1.4.0],[linbox-use@googlegroups.com],[linbox],
		[http://www.linalg.org/])
AC_CONFIG_MACRO_DIR([macros])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([1.8 gnu no-dependencies -Wall -Wno-portability])
AC_CONFIG_HEADERS([config.h])
AX_PREFIX_CONFIG_H(linbox/linbox-config.h, __LINBOX)
AC_PATH_PROG(RM, rm, $FALSE)
RM="$RM -f"


AM_MAINTAINER_MODE
AM_DEP_TRACK
AM_OUTPUT_DEPENDENCY_COMMANDS

AM_ACLOCAL_INCLUDE(macros)

# work around to fix the backward compatibility issue of automake 1.10 with 1.9 (pb with MKDIR_P)
AC_SUBST([MKDIR_P])

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
AX_CXX_COMPILE_STDCXX_11([noext],[optional])

STDFLAG=${CXXFLAGS}
CXXFLAGS=""


AC_COMPILER_NAME

AC_SUBST([DEFAULT_CFLAGS])
AC_SUBST([DEBUG_CFLAGS])
AC_SUBST([TESTS_CFLAGS])


TESTS_CFLAGS="-O0"
DEBUG_CFLAGS="-g"

DEFAULT_CFLAGS=""
# WARN_CFLAGS="-Wall"

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
		WARN_CFLAGS="${WARN_CFLAGS} -Wall -Wno-unused-parameter -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith  -Wwrite-strings -Wno-long-long"
		WARN_CFLAGS="${WARN_CFLAGS} -Wextra -ansi"
			;;
		xeko)
			WARN_CFLAGS="${WARN_CFLAGS} -Wno-unused-parameter"
			;;
		xgcc)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -ansi -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long"
			fi
			;;
		xgcc48)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -ansi -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long -fsanitize=address"
			fi
			;;

		xclang)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long"
			fi
			;;
		xclang31)
			WARN_CFLAGS="${WARN_CFLAGS} -Wextra -Wno-unused-parameter"
			if test "x${WARN}" = "xfull" ; then
				WARN_CFLAGS="${WARN_CFLAGS} -Wuninitialized -Wconversion -Wcast-qual -pedantic -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wno-long-long -fsanitize=address  -g"
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

# newer libtool...
LT_PREREQ([2.2])
LT_INIT



echo "-----------------------------------------------"


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

# Create some useful data types of fixed, known lengths

# We hereby assume that a character is always one byte
# LINBOX_INT8="char";

# case $ac_cv_sizeof_char in
  # 1)
     # TWO_BYTES=2;
     # FOUR_BYTES=4;
     # EIGHT_BYTES=8;
     # ;;
  # 8)
     # TWO_BYTES=16;
     # FOUR_BYTES=32;
     # EIGHT_BYTES=64;
# esac

# case $TWO_BYTES in
  # $ac_cv_sizeof_short)
     # LINBOX_INT16="short";
     # ;;
  # $ac_cv_sizeof_int)
     # LINBOX_INT16="int";
     # ;;
# esac

# case $FOUR_BYTES in
  # $ac_cv_sizeof_short)
     # LINBOX_INT32="short";
     # ;;
  # $ac_cv_sizeof_int)
     # LINBOX_INT32="int";
     # ;;
  # $ac_cv_sizeof_long)
     # LINBOX_INT32="long";
     # ;;
# esac

# case $EIGHT_BYTES in
  # $ac_cv_sizeof_short)
     # LINBOX_INT64="short";
     # ;;
  # $ac_cv_sizeof_int)
     # LINBOX_INT64="int";
     # ;;
  # $ac_cv_sizeof_long)
     # LINBOX_INT64="long";
     # ;;
  # $ac_cv_sizeof_long_long)
     # LINBOX_INT64="long long";
     # ;;
  # $ac_cv_sizeof___int64)
     # LINBOX_INT64="__int64";
     # ;;
# esac

# AC_DEFINE_UNQUOTED(INT8, $LINBOX_INT8, Canonical 8-bit data type)
# AC_DEFINE_UNQUOTED(INT16, $LINBOX_INT16, Canonical 16-bit data type)
# AC_DEFINE_UNQUOTED(INT32, $LINBOX_INT32, Canonical 32-bit data type)
# AC_DEFINE_UNQUOTED(INT64, $LINBOX_INT64, Canonical 64-bit data type)

echo "-----------------------------------------------"
# Feature checks
LB_MISC
LB_DRIVER


echo "-----------------------------------------------"

LB_CHECK_GMP(,,[
echo '*******************************************************************************'
echo ' ERROR: GMP not found!'
echo
echo ' GMP library compiled with --enable-cxx is required for this library to compile.'
echo ' Please make sure GMP is installed and specify its location with the option'
echo ' --with-gmp=<prefix> when running configure.'
echo '*******************************************************************************'
exit 1
])

LB_CHECK_IML
LB_CHECK_M4RI
LB_CHECK_M4RIE
LB_CHECK_MPFR
# we nedd mpfr next :
LB_CHECK_FPLLL
LB_CHECK_FLINT

LB_CHECK_NTL
LB_CHECK_GIVARO(,,[
echo '*******************************************************************************'
echo ' ERROR: GIVARO not found!'
echo
echo ' GIVARO library is required for this library to compile.'
echo ' Please make sure GIVARO is installed and specify its location with the'
echo ' option --with-givaro=<prefix> when running configure.'
echo '*******************************************************************************'
exit 1
])


LB_CHECK_SACLIB
LB_CHECK_LIDIA
LB_CHECK_MAPLE
LB_CHECK_SAGE
#LB_CHECK_ATLAS

# LB_CHECK_BLAS(,,[
# echo ''
# echo '*******************************************************************************'
# echo ' ERROR: BLAS not found!'
# echo
# echo ' BLAS routines are required for this library to compile. Please'
# echo ' make sure BLAS are installed and specify its location with the option'
# echo ' --with-blas=<lib> when running configure.'
# echo '*******************************************************************************'
# exit 1
# ])


LB_CHECK_FFLAS_FFPACK(,,[
echo ''
echo '*******************************************************************************'
echo ' ERROR: Fflas-Ffpack not found!'
echo
echo ' Fflas-Ffpack routines  are required for this library to compile. Please'
echo ' make sure they are installed and specify its location with the option'
echo ' --with-fflas-ffpack=<lib> when running configure.'
echo '*******************************************************************************'
exit 1
])

LB_CHECK_LAPACK

# AC_MSG_RESULT($BLAS_LIBS)
# AC_MSG_RESULT($BLAS_CFLAGS)


LB_CHECK_EXPAT

LB_BENCH
# LB_OPT
LB_DOC

LB_CHECK_OCL


if test ! -d ./benchmarks/data  ; then
	echo "Creating data dir in benchmark" ;
	mkdir ./benchmarks/data  ;
fi

# this is bogus so far as I can tell.  -bds
#CXXFLAGS="${GMP_CFLAGS} ${NTL_CFLAGS} ${GIVARO_CFLAGS} ${SACLIB_CFLAGS} ${CXXFLAGS}"
#AC_SUBST(CXXFLAGS)

CXXFLAGS="${GMP_CFLAGS} ${CXXFLAGS} ${STDFLAG}"

AC_SUBST(CXXFLAGS)

echo "-----------------------------------------------"

AC_CONFIG_FILES([
Makefile
linbox-config
examples/Makefile
doc/Makefile
linbox/Makefile
linbox/algorithms/Makefile
linbox/algorithms/gauss/Makefile
linbox/algorithms/IML/Makefile
linbox/algorithms/matrix-blas3/Makefile
linbox/algorithms/opencl-kernels/Makefile
linbox/blackbox/Makefile
linbox/element/Makefile
linbox/field/Makefile
linbox/field/NTL/Makefile
linbox/field/Givaro/Makefile
linbox/field/Modular/Makefile
linbox/matrix/Makefile
linbox/matrix/SparseMatrix/Makefile
linbox/randiter/Makefile
linbox/ring/Makefile
linbox/solutions/Makefile
linbox/switch/Makefile
linbox/util/Makefile
linbox/util/formats/Makefile
linbox/vector/Makefile
tests/Makefile
tests/data/Makefile
tests/matrix/Makefile
interfaces/Makefile
interfaces/driver/Makefile
interfaces/maple/Makefile
interfaces/kaapi/Makefile
interfaces/sage/Makefile
macros/Makefile
benchmarks/Makefile
])

# examples/Readme-make
# examples/field
# examples/solver
# examples/blackbox
# examples/data

AC_OUTPUT
