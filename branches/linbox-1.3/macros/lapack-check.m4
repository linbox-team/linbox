dnl Check for LAPACK in fflas-ffpack
dnl
dnl Copyright (c) 2011 the LinBox group
dnl Written by BB <bboyer@imag.fr>
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

dnl LB_CHECK_LAPACK

AC_DEFUN([LB_CHECK_LAPACK],
[

AC_MSG_CHECKING(for LAPACK in fflas-ffpack)

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${BACKUP_CXXFLAGS} ${FFLAS_FFPACK_CFLAGS} ${BLAS_CFLAGS}"
LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

AC_TRY_RUN(dnl ICC ?
[   #include "fflas-ffpack/fflas-ffpack-config.h"
	#ifdef __FFLASFFPACK_HAVE_LAPACK
	   int main() { return 0 ; }
   #else
   a pas lapack
   #endif
],dnl
[dnl
   AC_MSG_RESULT(ok)
   AC_DEFINE(HAVE_LAPACK,1,[Define if LAPACK is available])
   AM_CONDITIONAL(LINBOX_HAVE_LAPACK, true)
],dnl
[
	AC_MSG_RESULT(no)
	AM_CONDITIONAL(LINBOX_HAVE_LAPACK, false)
],dnl
[
   AC_MSG_RESULT(unknown)
   AM_CONDITIONAL(LINBOX_HAVE_LAPACK, false)
])

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}


])



