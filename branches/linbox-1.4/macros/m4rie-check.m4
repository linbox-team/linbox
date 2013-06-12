dnl Check for M4RIE
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
dnl Boyer Brice 18/04/12
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor

dnl LB_CHECK_M4RIE ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for M4RIE and define M4RIE_CFLAGS and M4RIE_LIBS

AC_DEFUN([LB_CHECK_M4RIE],
[

AC_ARG_WITH(m4rie,
[AC_HELP_STRING([--with-m4rie=<path>|yes], [Use M4RIE library. This library is (not yet) mandatory for
    LinBox compilation. If argument is yes or <empty> or <bad>
    that means the library is reachable with the standard
    search path (/usr or /usr/local). Otherwise you give
    the <path> to the directory which contains the
	library.  ])
])

AS_IF([test "$withval" = yes ],
	[ M4RIE_HOME_PATH="${DEFAULT_CHECKING_PATH}" ],
	[ test "$withval" != no ],
	[ M4RIE_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}" ],
	[ M4RIE_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl  min_m4rie_version=ifelse([$1], ,1.0.3,$1)


dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for M4RIE)

for M4RIE_HOME in ${M4RIE_HOME_PATH}
  do
    AS_IF([test -r "$M4RIE_HOME/include/m4rie/m4rie.h"],[

       AS_IF([ test "x$M4RIE_HOME" != "x/usr" -a "x$M4RIE_HOME" != "x/usr/local"],[
           M4RIE_CFLAGS="-I${M4RIE_HOME}/include"
           M4RIE_LIBS="-L${M4RIE_HOME}/lib -lm4rie"], [
		   M4RIE_CFLAGS=
		   M4RIE_LIBS="-lm4rie" ])

       CXXFLAGS="${BACKUP_CXXFLAGS} ${M4RIE_CFLAGS} ${GMP_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${M4RIE_LIBS} ${GMP_LIBS}"

       AC_TRY_LINK(
       [#include <m4rie/m4rie.h> ],
       [gf2e a;],
       [
	   AC_TRY_RUN(
	   [
	   int main () { return 0; /* not possible to check version */ }
	   ],[
	   m4rie_found="yes"
	   break
	   ],[
	   m4rie_problem="$problem $M4RIE_HOME"
	   unset M4RIE_CFLAGS
	   unset M4RIE_LIBS
	   ],[
	   m4rie_found="yes"
	   m4rie_cross="yes"
	   break
	   ])
	   ],
       [
       m4rie_found="no"
       m4rie_checked="$checked $M4RIE_HOME"
       unset M4RIE_CFLAGS
       unset M4RIE_LIBS
       ])
	   dnl  AC_MSG_RESULT(found in $m4rie_checked ? $m4rie_found)
    ],[
       m4rie_found="no"
	   dnl  AC_MSG_RESULT(not found at all $M4RIE_HOME : $m4rie_found)
    ])
done

AS_IF([test "x$m4rie_found" = "xyes"],[
		AC_SUBST(M4RIE_CFLAGS)
		AC_SUBST(M4RIE_LIBS)
		AC_DEFINE(HAVE_M4RIE,1,[Define if M4RIE is installed])
		HAVE_M4RIE=yes
		AS_IF([test "x$m4rie_cross" != "xyes"],
			[AC_MSG_RESULT(found)],
			[AC_MSG_RESULT(unknown)
			echo "WARNING: You appear to be cross compiling, so there is no way to determine"
			echo "whether your M4RIE version is new enough. I am assuming it is."
			])
		ifelse([$2], , :, [$2])],
		[ test -n "$m4rie_problem" ],
		[ AC_MSG_RESULT(problem)
		echo "Sorry, your M4RIE version is too old. Disabling."
		ifelse([$3], , :, [$3]) ],
		[ test "x$m4rie_found" = "xno" ],
		[ AC_MSG_RESULT(not found)
		ifelse([$3], , :, [$3])])

AM_CONDITIONAL(LINBOX_HAVE_M4RIE, test "x$HAVE_M4RIE" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
