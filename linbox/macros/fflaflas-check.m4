# Check for Fflas-Ffpack
# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence
# Boyer Brice 19/04/11
# Bradford Hovinen, 2001-06-13
# Modified by Pascal Giorgi, 2003-12-03
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_FFLAFLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for FFLAFLAS and define FFLAFLAS_CFLAGS and FFLAFLAS_LIBS

AC_DEFUN([LB_CHECK_FFLAFLAS],
[

AC_ARG_WITH(fflas-ffpack,
	[  --with-fflas-ffpack=<path>|yes Use Fflas-Ffpack library. This library is mandatory for
		LinBox compilation. If argument is yes or <empty> or <bad> :)
		that means the library is reachable with the standard
		search path (/usr or /usr/local). Otherwise you give
		the <path> to the directory which contains the
		library.
	],
    [if test "$withval" = yes ; then
        FFLAFLAS_HOME_PATH="${DEFAULT_CHECKING_PATH}"
        elif test "$withval" != no ; then
        FFLAFLAS_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
        fi],
    [FFLAFLAS_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl  min_iml_version=ifelse([$1], ,1.0.3,$1)


dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for FFLAFLAS)

for FFLAFLAS_HOME in ${FFLAFLAS_HOME_PATH}
  do
    if test -r "$FFLAFLAS_HOME/include/fflas-ffpack/fflas-ffpack.h"; then

		BLAS_LIBS=`$FFLAFLAS_HOME/bin/fflasffpack-config --blas-libs`
		dnl  BLAS_CFLAGS=`$FFLAFLAS_HOME/bin/fflasffpack-config --cflags`


       if test "x$FFLAFLAS_HOME" != "x/usr" -a "x$FFLAFLAS_HOME" != "x/usr/local"; then
           FFLAFLAS_CFLAGS="-I${FFLAFLAS_HOME}/include"
       else
           FFLAFLAS_CFLAGS=
       fi

       CXXFLAGS="${BACKUP_CXXFLAGS} ${FFLAFLAS_CFLAGS} ${BLAS_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

       AC_TRY_LINK(
       [#include "fflas-ffpack/fflas-ffpack.h"],
       [FFLAS::FFLAS_TRANSPOSE a;],
       [
	   fflaflas_found="yes"
	   FFLAFLAS_LOC="$FFLAFLAS_HOME"
	   ],
       [
       fflaflas_found="no"
       fflaflas_checked="$checked $FFLAFLAS_HOME"
       unset FFLAFLAS_CFLAGS
	   unset FFLAFLAS_LOC
	   unset BLAS_LIBS
       ])
	   dnl  AC_MSG_RESULT(found in $fflaflas_checked ? $fflaflas_found)
    else
       fflasflas_found="no"
	   dnl  AC_MSG_RESULT(not found at all $FFLAFLAS_HOME : $fflaflas_found)
    fi
done

if test "x$fflaflas_found" = "xyes" ; then
    AC_SUBST(FFLAFLAS_CFLAGS)
    AC_SUBST(FFLAFLAS_LIBS)
	AC_SUBST(FFLAFLAS_LOC)
	AC_SUBST(BLAS_LIBS)
    AC_DEFINE(HAVE_FFLAFLAS,1,[Define if FFLAFLAS is installed])
    HAVE_FFLAFLAS=yes
    if test "x$fflasflas_cross" != "xyes"; then
        AC_MSG_RESULT(found)
    else
        AC_MSG_RESULT(unknown)
        echo "WARNING: You appear to be cross compiling, so there is no way to determine"
        echo "whether your FFLAFLAS version is new enough. I am assuming it is."
    fi
    ifelse([$2], , :, [$2])
elif test -n "$fflasflas_problem"; then
    AC_MSG_RESULT(problem)
    echo "Sorry, your FFLAFLAS version is too old. Disabling."
    ifelse([$3], , :, [$3])
elif test "x$fflasflas_found" = "xno" ; then
    AC_MSG_RESULT(not found)
    ifelse([$3], , :, [$3])
fi

AM_CONDITIONAL(LINBOX_HAVE_FFLAFLAS, test "x$HAVE_FFLAFLAS" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
