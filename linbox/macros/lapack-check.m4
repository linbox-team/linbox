dnl Check for ATLAS
dnl Copyright (c) 2011 the LinBox group
dnl Written by BB <bboyer@imag.fr>
dnl This file is part of LinBox
dnl see COPYING for licence

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
	#ifdef __FFLAS_FFPACK_HAVE_LAPACK
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



