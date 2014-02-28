dnl Check for FLINT
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
dnl Boyer Brice 8/6/12
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor

dnl LB_CHECK_FLINT ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for FLINT and define FLINT_CFLAGS and FLINT_LIBS

AC_DEFUN([LB_CHECK_FLINT],
[

AC_ARG_WITH(flint,
[AC_HELP_STRING([--with-flint=<path>|yes], [Use FLINT library. This library is (not yet) mandatory for
    LinBox compilation. If argument is yes or <empty> or <bad>
    that means the library is reachable with the standard
    search path (/usr or /usr/local). Otherwise you give
    the <path> to the directory which contains the
	library.  ])
])

AS_IF([test "$withval" = yes ],
	[ FLINT_HOME_PATH="${DEFAULT_CHECKING_PATH}" ],
	[ test "$withval" != no ],
	[ FLINT_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}" ],
	[ FLINT_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl  min_flint_version=ifelse([$1], ,1.0.3,$1)


dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for FLINT)

for FLINT_HOME in ${FLINT_HOME_PATH}
  do
    AS_IF([test -r "$FLINT_HOME/include/flint/flint.h"],[

       AS_IF([ test "x$FLINT_HOME" != "x/usr" -a "x$FLINT_HOME" != "x/usr/local"],[
           FLINT_CFLAGS="-I${FLINT_HOME}/include"
           FLINT_LIBS="-L${FLINT_HOME}/lib -lflint"], [
		   FLINT_CFLAGS=
		   FLINT_LIBS="-lflint" ])

       CXXFLAGS="${BACKUP_CXXFLAGS} ${GMP_CFLAGS} ${MPFR_CFLAGS} ${FLINT_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${GMP_LIBS} ${MPFR_LIBS} ${FLINT_LIBS}"

       AC_TRY_LINK(
       [ //extern "C" {
      #define __GMP_BITS_PER_MP_LIMB GMP_LIMB_BITS
       #include <flint/flint.h>
       #include <flint/fmpz_mat.h>
       //}
       ],
       [fmpz_mat_t a;],
       [
	   AC_TRY_RUN(
	   [
	   #include "flint/flint.h"
	   int main () {
	   if ( (__FLINT_VERSION < 2) || ( (__FLINT_VERSION ==2) && (__FLINT_VERSION_MINOR < 4) ) )
	     return 1 ;
	   else
	     return 0;
	   }
	   ],[
	   flint_found="yes"
	   break
	   ],[
	   flint_problem="$problem $FLINT_HOME"
	   unset FLINT_CFLAGS
	   unset FLINT_LIBS
	   ],[
	   flint_found="yes"
	   flint_cross="yes"
	   break
	   ])
	   ],
       [
       flint_found="no"
       flint_checked="$checked $FLINT_HOME"
       unset FLINT_CFLAGS
       unset FLINT_LIBS
       ])
	   dnl  AC_MSG_RESULT(found in $flint_checked ? $flint_found)
    ],[
       flint_found="no"
	   dnl  AC_MSG_RESULT(not found at all $FLINT_HOME : $flint_found)
    ])
done

AS_IF([test "x$flint_found" = "xyes"],[
		AC_SUBST(FLINT_CFLAGS)
		AC_SUBST(FLINT_LIBS)
		AC_DEFINE(HAVE_FLINT,1,[Define if FLINT is installed])
		HAVE_FLINT=yes
		AS_IF([test "x$flint_cross" != "xyes"],
			[AC_MSG_RESULT(found)],
			[AC_MSG_RESULT(unknown)
			echo "WARNING: You appear to be cross compiling, so there is no way to determine"
			echo "whether your FLINT version is new enough. I am assuming it is."
			])
		ifelse([$2], , :, [$2])],
		[ test -n "$flint_problem" ],
		[ AC_MSG_RESULT(problem)
		echo "Sorry, your FLINT version is too old. Disabling."
		ifelse([$3], , :, [$3]) ],
		[ test "x$flint_found" = "xno" ],
		[ AC_MSG_RESULT(>=2.4 not found )
		ifelse([$3], , :, [$3])])

AM_CONDITIONAL(LINBOX_HAVE_FLINT, test "x$HAVE_FLINT" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
