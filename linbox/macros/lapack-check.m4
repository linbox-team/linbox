# Check for ATLAS
# Copyright (c) 2011 the LinBox group
# Written by BB <bboyer@imag.fr>
# This file is part of LinBox
# see COPYING for licence

dnl LB_CHECK_LAPACK

AC_DEFUN([LB_CHECK_LAPACK],
[

AC_MSG_CHECKING(for LAPACK in fflas-ffpack)

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${BACKUP_CXXFLAGS} ${FFLAFLAS_CFLAGS} ${BLAS_CFLAGS}"
LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

AC_TRY_RUN(dnl ICC ?
[   #include "fflasffpack-config.h"
	#ifdef __FFLAFLAS_HAVE_LAPACK
	   int main() { return 0 ; }
   #else
   a pas lapack
   #endif
],dnl
[dnl
   AC_MSG_RESULT(ok)
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



