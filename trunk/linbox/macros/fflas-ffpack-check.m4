dnl Check for Fflas-Ffpack
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
dnl Boyer Brice 19/04/11
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor

dnl LB_CHECK_FFLAS_FFPACK ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for FFLAS_FFPACK and define FFLAS_FFPACK_CFLAGS and FFLAS_FFPACK_LIBS

AC_DEFUN([LB_CHECK_FFLAS_FFPACK],
[

AC_ARG_WITH(fflas-ffpack,
	[AC_HELP_STRING([--with-fflas-ffpack=<path>|yes], [Use Fflas-Ffpack library. This library is mandatory for
		LinBox compilation. If argument is yes or <empty> or <bad> :)
		that means the library is reachable with the standard
		search path (/usr or /usr/local). Otherwise you give
		the <path> to the directory which contains the
		library.
        Available at "http://linalg.org/projects/fflas-ffpack".
	])],
    [if test "$withval" = yes ; then
        FFLAS_FFPACK_HOME_PATH="${DEFAULT_CHECKING_PATH}"
        elif test "$withval" != no ; then
        FFLAS_FFPACK_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
        fi],
    [FFLAS_FFPACK_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl  min_iml_version=ifelse([$1], ,1.0.3,$1)

dnl -------------------- dnl
dnl FFLAS-FFPACK VERSION dnl
dnl -------------------- dnl

version_min=10700
version_max=10800


dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for FFLAS-FFPACK >= $version_min and < $version_max)

for FFLAS_FFPACK_HOME in ${FFLAS_FFPACK_HOME_PATH}
  do
    if test -r "$FFLAS_FFPACK_HOME/include/fflas-ffpack/fflas-ffpack.h" -a -x "$FFLAS_FFPACK_HOME/bin/fflas-ffpack-config"; then

		BLAS_LIBS=`$FFLAS_FFPACK_HOME/bin/fflas-ffpack-config --libs`
		BLAS_CFLAGS=`$FFLAS_FFPACK_HOME/bin/fflas-ffpack-config --cflags`


       if test "x$FFLAS_FFPACK_HOME" != "x/usr" -a "x$FFLAS_FFPACK_HOME" != "x/usr/local"; then
           FFLAS_FFPACK_CFLAGS="-I${FFLAS_FFPACK_HOME}/include"
       else
           FFLAS_FFPACK_CFLAGS=
       fi

       CXXFLAGS="${BACKUP_CXXFLAGS} ${FFLAS_FFPACK_CFLAGS} ${BLAS_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${BLAS_LIBS}"

       AC_TRY_LINK(
       [#include "fflas-ffpack/fflas-ffpack.h"],
       [FFLAS::FFLAS_TRANSPOSE a;],
       [

	FF_VER=`$FFLAS_FFPACK_HOME/bin/fflas-ffpack-config --decimal-version`
	AS_IF([ test $FF_VER -ge $version_min -a $FF_VER -lt $version_max ],
		[
		ffflasffpack_found="yes"
		FFLAS_FFPACK_LOC="$FFLAS_FFPACK_HOME"
		break
		],
		[
		ffflasffpack_found="no"
		]
		)
	],
       [
       ffflasffpack_found="no"
       ffflasffpack_checked="$checked $FFLAS_FFPACK_HOME"
       unset FFLAS_FFPACK_CFLAGS
	   unset FFLAS_FFPACK_LOC
	   unset BLAS_LIBS
	   unset BLAS_CFLAGS
       ])
	   dnl  AC_MSG_RESULT(found in $ffflasffpack_checked ? $ffflasffpack_found)
    else
       fflasflas_found="no"
	   dnl  AC_MSG_RESULT(not found at all $FFLAS_FFPACK_HOME : $ffflasffpack_found)
    fi
done

if test "x$ffflasffpack_found" = "xyes" ; then
    AC_SUBST(FFLAS_FFPACK_CFLAGS)
    AC_SUBST(FFLAS_FFPACK_LIBS)
	AC_SUBST(FFLAS_FFPACK_LOC)
	AC_SUBST(BLAS_LIBS)
	AC_SUBST(BLAS_CFLAGS)
    AC_DEFINE(HAVE_FFLAS_FFPACK,1,[Define if FFLAS-FFPACK is installed])
	FF_VER=`$FFLAS_FFPACK_LOC/bin/fflas-ffpack-config --decimal-version`
	AC_DEFINE_UNQUOTED(FFLAS_FFPACK_VERSION, $FF_VER ,[what version of FFLAS-FFPACK is installed])
	HAVE_FFLAS_FFPACK=yes
    if test "x$fflasflas_cross" != "xyes"; then
        AC_MSG_RESULT(found)
    else
        AC_MSG_RESULT(unknown)
        echo "WARNING: You appear to be cross compiling, so there is no way to determine"
        echo "whether your FFLAS-FFPACK version is new enough. I am assuming it is."
    fi
    ifelse([$2], , :, [$2])
elif test -n "$fflasflas_problem"; then
    AC_MSG_RESULT(problem)
    echo "Sorry, your FFLAS-FFPACK version is too old. Disabling."
    ifelse([$3], , :, [$3])
elif test "x$fflasflas_found" = "xno" ; then
    AC_MSG_RESULT(not found)
    ifelse([$3], , :, [$3])
fi

AM_CONDITIONAL(LINBOX_HAVE_FFLAS_FFPACK, test "x$HAVE_FFLAS_FFPACK" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
