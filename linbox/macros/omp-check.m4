dnl turn on OPENMP
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

dnl LB_CHECK_OMP()
dnl
dnl turn on OpenMP

AC_DEFUN([LB_OMP],[
	AC_MSG_CHECKING(enabling OpenMP)

AC_ARG_WITH(
	[openmp],
	[AS_HELP_STRING([--with-openmp],[Turn on OpenMP based parallel code])],
	[with_openmp=yes],
	[with_openmp=no])

AS_IF([ test "$with_openmp" = "yes" ],
	[AC_MSG_RESULT(yes)
	CXXFLAGS="-fopenmp -D__FFLASFFPACK_USE_OPENMP ${CXXFLAGS}"
	],[
	AC_MSG_RESULT(no)
	])
])

