dnl Check for GIVARO
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
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor

dnl LB_CHECK_GIVARO ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for Givaro and define GIVARO_CFLAGS and GIVARO_LIBS

AC_DEFUN([LB_CHECK_GIVARO],
[

AC_ARG_WITH(givaro,
[AC_HELP_STRING([--with-givaro=<path>|yes], [Use Givaro library. This library is mandatory for
                           LinBox compilation. If argument is yes or <empty>
			   that means the library is reachable with the standard
			   search path (/usr or /usr/local). Otherwise you give
			   the <path> to the directory which contains the
			   library.
])],
	     [if test "$withval" = yes ; then
			GIVARO_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			GIVARO_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
	     fi],
	     [GIVARO_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl -------------- dnl
dnl GIVARO VERSION dnl
dnl -------------- dnl

version_min=30700
version_max=30800

dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for GIVARO >= $version_min and < $version_max)

for GIVARO_HOME in ${GIVARO_HOME_PATH}
 do
if test -r "$GIVARO_HOME/include/givaro/givconfig.h"; then

	if test "x$GIVARO_HOME" != "x/usr" -a "x$GIVARO_HOME" != "x/usr/local"; then
		GIVARO_CFLAGS="-I${GIVARO_HOME}/include"
		GIVARO_LIBS="-L${GIVARO_HOME}/lib -lgivaro"
	else
		GIVARO_CFLAGS=
		GIVARO_LIBS="-lgivaro"
	fi
	CXXFLAGS="${BACKUP_CXXFLAGS} ${GIVARO_CFLAGS} ${GMP_CFLAGS}"
	LIBS="${BACKUP_LIBS} ${GIVARO_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK(
	[#include <givaro/givinteger.h>],
	[Givaro::Integer a;],
	[
	AC_TRY_RUN(
	[#include <givaro/givconfig.h>
	 int main () { if (GIVARO_VERSION < $version_min || GIVARO_VERSION >= $version_max || GIVARO_VERSION>0x030000) return -1; else return 0; /* old version of Givaro are defined as hexa 0x03yyzz*/ }
	],[
	givaro_found="yes"
	break
	],[
	givaro_problem="$problem $GIVARO_HOME"
	unset GIVARO_CFLAGS
	unset GIVARO_LIBS
	],[
	givaro_found="yes"
	givaro_cross="yes"

	break
	])
	],
	[
	givaro_found="no"
	givaro_checked="$checked $GIVARO_HOME"
	unset GIVARO_CFLAGS
	unset GIVARO_LIBS

	])
else
	givaro_found="no"
fi
done

if test "x$givaro_found" = "xyes" ; then
	AC_SUBST(GIVARO_CFLAGS)
	AC_SUBST(GIVARO_LIBS)
	dnl  echo $GIVARO_CFLAGS $GIVARO_LIBS
	AC_DEFINE(HAVE_GIVARO,1,[Define if GIVARO is installed])
	HAVE_GIVARO=yes
	if test "x$givaro_cross" != "xyes"; then
		AC_MSG_RESULT(found)
	else
		AC_MSG_RESULT(unknown)
		echo "WARNING: You appear to be cross compiling, so there is no way to determine"
		echo "whether your GIVARO version is new enough. I am assuming it is."
	fi
	ifelse([$2], , :, [$2])
elif test -n "$givaro_problem"; then
	AC_MSG_RESULT(problem)
	echo "Sorry, your GIVARO version is too old. Disabling."
	ifelse([$3], , :, [$3])
elif test "x$givaro_found" = "xno" ; then
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi

AM_CONDITIONAL(LINBOX_HAVE_GIVARO, test "x$HAVE_GIVARO" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
