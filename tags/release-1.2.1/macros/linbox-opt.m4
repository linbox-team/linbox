# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence



AC_DEFUN([LB_OPT],
[
AC_MSG_CHECKING([whether to use run time optimization])

AC_ARG_ENABLE(optimization,
[AC_HELP_STRING([--enable-optimization], [Enable run time optimization in LinBox code (only Strassen matrix threshold for now)])],
[


BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS=${FFLAFLAS_CFLAGS}

AC_TRY_RUN([ #include "fflasffpack-config.h"
   int main() {
#ifdef __FFLAFLAS_STRASSEN_OPTIMIZATION
return 0;
#else
pas bon !
#endif
} ],[dnl OK
	strassen_opti="yes"
	AC_DEFINE(STRASSEN_OPTIMIZATION,,[Define if optimized  threshold for Strassen-Winograd matrix multiplication is available])
    WINO="`grep "define.*__FFLAFLAS_WINOTHRESHOLD" ${FFLAFLAS_LOC}/include/fflasffpack-config.h  | awk '{print $NF}'`"
	AC_MSG_RESULT(ok : $WINO)
	AC_DEFINE_UNQUOTED(WINOTHRESHOLD, $WINO, [optimized threshold for switching to strassen matrix multiplication])
	],[ dnl NO
	AC_MSG_RESULT(not enabled. Please optimise Strassen threshold in Fflas-Ffpack)
	strassen_opti="no"
	],[ dnl CROSS
	AC_MSG_RESULT(cross compilation)
	strassen_opti="no"
	])

],
[AC_MSG_RESULT(no)])
])

CXXFLAGS=${BACKUP_CXXFLAGS}