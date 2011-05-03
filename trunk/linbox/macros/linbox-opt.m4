# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence



AC_DEFUN([LB_OPT],
[
AC_MSG_CHECKING([whether to use run time optimization])

AC_ARG_ENABLE(optimization,
[AC_HELP_STRING([--enable-optimization], [Enable run time optimization in LinBox code (only Strassen matrix threshold for now)])],
[
AC_MSG_RESULT(yes)


BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

if test "x$HAVE_BLAS" = "xyes" ;then
AC_MSG_CHECKING([best threshold for Strassen-Winograd matrix multiplication])


CXXFLAGS="${BACKUP_CXXFLAGS} -I`pwd` -I`pwd`/linbox ${BLAS_CFLAGS} ${FFLAFLAS_CFLAGS} ${GMP_CFLAGS}  ${GIVARO_CFLAGS} ${CBLAS_FLAG}"
LIBS="${BACKUP_LIBS} ${BLAS_LIBS} ${GIVARO_LIBS} ${GMP_LIBS} "


echo   " #define __LINBOX_INT8  $LINBOX_INT8
	 #define __LINBOX_INT16 $LINBOX_INT16
	 #define __LINBOX_INT32 $LINBOX_INT32
	 #define __LINBOX_INT64 $LINBOX_INT64
" > linbox/linbox-config.h


AC_TRY_RUN([ #include "fflas-ffpack/fflasffpack-config.h"
   int main() {
#ifdef __FFLAFLAS_STRASSEN_OPTIMIZATION
return 0;
#else
pas bon !
#endif
} ],[
	],[
	strassen_opti="yes"
	 AC_DEFINE(STRASSEN_OPTIMIZATION,,[Define if optimized  threshold for Strassen-Winograd matrix multiplication is available])
     WINO="`grep "define.*__FFLAFLAS_WINOTHRESHOLD" ${FFLAFLAS_LOC}/include/fflasffpack-config.h  | awk '{print $NF}'`"
	 AC_MSG_RESULT(ok : $WINO)
	 AC_DEFINE_UNQUOTED(WINOTHRESHOLD, $WINO, [optimized threshold for switching to strassen matrix multiplication])
	fi
	],[
	AC_MSG_RESULT(not enabled. Please optimise Strassen threshold in Fflas-Ffpack)
	strassen_opti="no"
	],[
	AC_MSG_RESULT(cross compilation)
	strassen_opti="no"
	])

])
])
