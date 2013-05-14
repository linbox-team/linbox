dnl Check for FPLLL
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
dnl Boyer Brice 22/10/11
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor

dnl LB_CHECK_FPLLL ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for FPLLL and define FPLLL_CFLAGS and FPLLL_LIBS

AC_DEFUN([LB_CHECK_FPLLL],
[

AC_ARG_WITH(fplll,
[AC_HELP_STRING([--with-fplll=<path>|yes], [Use FPLLL library. This library is (not yet) mandatory for
    LinBox compilation. If argument is yes or <empty>
    that means the library is reachable with the standard
    search path (/usr or /usr/local). Otherwise you give
    the <path> to the directory which contains the
    library.
])],
    [if test "$withval" = yes ; then
        FPLLL_HOME_PATH="${DEFAULT_CHECKING_PATH}"
        elif test "$withval" != no ; then
        FPLLL_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
        fi],
    [FPLLL_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl  min_fplll_version=ifelse([$1], ,1.0.3,$1)


dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

version_min=4.0.1
AC_MSG_CHECKING(for FPLLL >= $version_min)


for FPLLL_HOME in ${FPLLL_HOME_PATH}
  do
    if test -r "$FPLLL_HOME/include/fplll/fplll.h"; then

       AS_IF([ test "x$FPLLL_HOME" != "x/usr" -a "x$FPLLL_HOME" != "x/usr/local"], [
           FPLLL_CFLAGS="-I${FPLLL_HOME}/include"
           FPLLL_LIBS="-L${FPLLL_HOME}/lib -lfplll"
       ],[
           FPLLL_CFLAGS=
           FPLLL_LIBS="-lfplll"
       ])

       CXXFLAGS="${BACKUP_CXXFLAGS} ${MPFR_CFLAGS} ${FPLLL_CFLAGS} ${GMP_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${GMP_LIBS} ${MPFR_LIBS} ${FPLLL_LIBS} "

       AC_TRY_LINK(
       [
	   #include <fplll/fplll.h>
	   ],
       [enum fplll::MatPrintMode a ;],
       [
	   AC_TRY_RUN(
	   [
	   int main () { return 0; /* not important to check for  version */ }
	   ],[
	   fplll_found="yes"
	   break
	   ],[
	   fplll_problem="$problem $FPLLL_HOME"
	   unset FPLLL_CFLAGS
	   unset FPLLL_LIBS
	   ],[
	   fplll_found="yes"
	   fplll_cross="yes"
	   break
	   ])
	   ],
       [
       fplll_found="no"
       fplll_checked="$checked $FPLLL_HOME"
       unset FPLLL_CFLAGS
       unset FPLLL_LIBS
       ])
	   dnl  AC_MSG_RESULT(found in $fplll_checked ? $fplll_found)
    else
       fplll_found="no"
	   dnl  AC_MSG_RESULT(not found at all $FPLLL_HOME : $fplll_found)
    fi
done

if test "x$fplll_found" = "xyes" ; then
    AC_SUBST(FPLLL_CFLAGS)
    AC_SUBST(FPLLL_LIBS)
    AC_DEFINE(HAVE_FPLLL,1,[Define if FPLLL is installed])
    HAVE_FPLLL=yes
    AS_IF([ test "x$fplll_cross" != "xyes" ],[
        AC_MSG_RESULT(found)
    ],[
        AC_MSG_RESULT(unknown)
        echo "WARNING: You appear to be cross compiling, so there is no way to determine"
        echo "whether your FPLLL version is new enough. I am assuming it is."
    ])
    ifelse([$2], , :, [$2])
elif test -n "$fplll_problem"; then
    AC_MSG_RESULT(problem)
    echo "Sorry, your FPLLL version is too old. Disabling."
    ifelse([$3], , :, [$3])
elif test "x$fplll_found" = "xno" ; then
    AC_MSG_RESULT(not found)
    ifelse([$3], , :, [$3])
fi

AM_CONDITIONAL(LINBOX_HAVE_FPLLL, test "x$HAVE_FPLLL" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
