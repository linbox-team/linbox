dnl Copyright (c) the LinBox group
dnl This file is part of LinBox

 dnl ========LICENCE========
 dnl This file is part of the library LinBox.
 dnl
 dnl LinBox is free software: you can redistribute it and/or modify
 dnl it under the terms of the  GNU Lesser General Public
 dnl License as published by the Free Software Foundation; either
 dnl version 2.1 of the License, or (at your option) any later version.
 dnl
 dnl This library is distributed in the hope that it will be useful,
 dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
 dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 dnl Lesser General Public License for more details.
 dnl
 dnl You should have received a copy of the GNU Lesser General Public
 dnl License along with this library; if not, write to the Free Software
 dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 dnl ========LICENCE========
 dnl



AC_DEFUN([LB_OPT],
[
AC_MSG_CHECKING([whether to use run time optimization])

AC_ARG_ENABLE(optimization,
[AC_HELP_STRING([--enable-optimization], [Enable run time optimization in LinBox code (only Strassen matrix threshold for now)])],
[
	dnl AC_DEFINE_UNQUOTED(WINOTHRESHOLD, __FFLASFFPACK_WINOTHRESHOLD) ?

	WINO="`grep "define.*__FFLASFFPACK_WINOTHRESHOLD" ${FFLAS_FFPACK_LOC}/include/fflas-ffpack/fflas-ffpack-optimise.h  | awk '{print $NF}'`"
	AS_IF([ test -z "{WINO}"],[
		AC_MSG_RESULT("fflas-ffpack was not optimised. Defaulting")
		WINO="`grep "define.*__FFLASFFPACK_WINOTHRESHOLD" ${FFLAS_FFPACK_LOC}/include/fflas-ffpack/fflas-ffpack-config.h  | awk '{print $NF}'`"
	],[
	AC_MSG_RESULT(OK)
	])
	AC_DEFINE_UNQUOTED(WINOTHRESHOLD, $WINO, [optimized threshold for switching to strassen matrix multiplication])

],
[AC_MSG_RESULT(no)])
])

