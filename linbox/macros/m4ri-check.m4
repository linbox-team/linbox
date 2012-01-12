dnl Check for M4RI
dnl Copyright (c) the LinBox group
dnl This file is part of LinBox
dnl see COPYING for licence
dnl Boyer Brice 22/06/11
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor
dnl This file is part of LinBox, see COPYING for licence information.

dnl LB_CHECK_M4RI ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for M4RI and define M4RI_CFLAGS and M4RI_LIBS

AC_DEFUN([LB_CHECK_M4RI],
[

AC_ARG_WITH(m4ri,
[AC_HELP_STRING([--with-m4ri=<path>|yes], [Use M4RI library. This library is (not yet) mandatory for
    LinBox compilation. If argument is yes or <empty> or <bad>
    that means the library is reachable with the standard
    search path (/usr or /usr/local). Otherwise you give
    the <path> to the directory which contains the
	library.  ])
])

AS_IF([test "$withval" = yes ],
	[ M4RI_HOME_PATH="${DEFAULT_CHECKING_PATH}" ],
	[ test "$withval" != no ],
	[ M4RI_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}" ],
	[ M4RI_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl  min_m4ri_version=ifelse([$1], ,1.0.3,$1)


dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for M4RI)

for M4RI_HOME in ${M4RI_HOME_PATH}
  do
    AS_IF([test -r "$M4RI_HOME/include/m4ri/m4ri.h"],[

       AS_IF([ test "x$M4RI_HOME" != "x/usr" -a "x$M4RI_HOME" != "x/usr/local"],[
           M4RI_CFLAGS="-I${M4RI_HOME}/include"
           M4RI_LIBS="-L${M4RI_HOME}/lib -lm4ri"], [
		   M4RI_CFLAGS=
		   M4RI_LIBS="-lm4ri" ])

       CXXFLAGS="${BACKUP_CXXFLAGS} ${M4RI_CFLAGS} ${GMP_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${M4RI_LIBS} ${GMP_LIBS}"

       AC_TRY_LINK(
       [#include <m4ri/m4ri.h> ],
       [mzd_t a;],
       [
	   AC_TRY_RUN(
	   [
	   int main () { return 0; /* not possible to check version */ }
	   ],[
	   m4ri_found="yes"
	   break
	   ],[
	   m4ri_problem="$problem $M4RI_HOME"
	   unset M4RI_CFLAGS
	   unset M4RI_LIBS
	   ],[
	   m4ri_found="yes"
	   m4ri_cross="yes"
	   break
	   ])
	   ],
       [
       m4ri_found="no"
       m4ri_checked="$checked $M4RI_HOME"
       unset M4RI_CFLAGS
       unset M4RI_LIBS
       ])
	   dnl  AC_MSG_RESULT(found in $m4ri_checked ? $m4ri_found)
    ],[
       m4ri_found="no"
	   dnl  AC_MSG_RESULT(not found at all $M4RI_HOME : $m4ri_found)
    ])
done

AS_IF([test "x$m4ri_found" = "xyes"],[
		AC_SUBST(M4RI_CFLAGS)
		AC_SUBST(M4RI_LIBS)
		AC_DEFINE(HAVE_M4RI,1,[Define if M4RI is installed])
		HAVE_M4RI=yes
		AS_IF([test "x$m4ri_cross" != "xyes"],
			[AC_MSG_RESULT(found)],
			[AC_MSG_RESULT(unknown)
			echo "WARNING: You appear to be cross compiling, so there is no way to determine"
			echo "whether your M4RI version is new enough. I am assuming it is."
			])
		ifelse([$2], , :, [$2])],
		[ test -n "$m4ri_problem" ],
		[ AC_MSG_RESULT(problem)
		echo "Sorry, your M4RI version is too old. Disabling."
		ifelse([$3], , :, [$3]) ],
		[ test "x$m4ri_found" = "xno" ],
		[ AC_MSG_RESULT(not found)
		ifelse([$3], , :, [$3])])

AM_CONDITIONAL(LINBOX_HAVE_M4RI, test "x$HAVE_M4RI" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
