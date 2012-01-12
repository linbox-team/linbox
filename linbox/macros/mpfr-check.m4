dnl Check for MPFR
dnl Copyright (c) the LinBox group
dnl This file is part of LinBox
dnl see COPYING for licence
dnl Boyer Brice 22/10/11
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor
dnl This file is part of LinBox, see COPYING for licence information.

dnl LB_CHECK_MPFR ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for MPFR and define MPFR_CFLAGS and MPFR_LIBS

AC_DEFUN([LB_CHECK_MPFR],
[

AC_ARG_WITH(mpfr,
[AC_HELP_STRING([--with-mpfr=<path>|yes], [Use MPFR library. This library is (not yet) mandatory for
    LinBox compilation. If argument is yes or <empty> or <bad> :)
    that means the library is reachable with the standard
    search path (/usr or /usr/local). Otherwise you give
    the <path> to the directory which contains the
    library.
])],
    [if test "$withval" = yes ; then
        MPFR_HOME_PATH="${DEFAULT_CHECKING_PATH}"
        elif test "$withval" != no ; then
        MPFR_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
        fi],
    [MPFR_HOME_PATH="${DEFAULT_CHECKING_PATH}"])

dnl  min_mpfr_version=ifelse([$1], ,1.0.3,$1)


dnl Check for existence
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for MPFR)

for MPFR_HOME in ${MPFR_HOME_PATH}
  do
    if test -r "$MPFR_HOME/include/mpfr.h"; then

       AS_IF([ test "x$MPFR_HOME" != "x/usr" -a "x$MPFR_HOME" != "x/usr/local"], [
           MPFR_CFLAGS="-I${MPFR_HOME}/include"
           MPFR_LIBS="-L${MPFR_HOME}/lib -lmpfr"
       ],[
           MPFR_CFLAGS=
           MPFR_LIBS="-lmpfr"
       ])

       CXXFLAGS="${BACKUP_CXXFLAGS} ${MPFR_CFLAGS} ${GMP_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${GMP_LIBS} ${MPFR_LIBS} "

       AC_TRY_LINK(
       [
	   #include <mpfr.h>
	   ],
       [mpfr_t a ;],
       [
	   AC_TRY_RUN(
	   [
	   int main () { return 0; /* not important to check for  version */ }
	   ],[
	   mpfr_found="yes"
	   break
	   ],[
	   mpfr_problem="$problem $MPFR_HOME"
	   unset MPFR_CFLAGS
	   unset MPFR_LIBS
	   ],[
	   mpfr_found="yes"
	   mpfr_cross="yes"
	   break
	   ])
	   ],
       [
       mpfr_found="no"
       mpfr_checked="$checked $MPFR_HOME"
       unset MPFR_CFLAGS
       unset MPFR_LIBS
       ])
	   dnl  AC_MSG_RESULT(found in $mpfr_checked ? $mpfr_found)
    else
       mpfr_found="no"
	   dnl  AC_MSG_RESULT(not found at all $MPFR_HOME : $mpfr_found)
    fi
done

if test "x$mpfr_found" = "xyes" ; then
    AC_SUBST(MPFR_CFLAGS)
    AC_SUBST(MPFR_LIBS)
    AC_DEFINE(HAVE_MPFR,1,[Define if MPFR is installed])
    HAVE_MPFR=yes
    AS_IF([ test "x$mpfr_cross" != "xyes" ],[
        AC_MSG_RESULT(found)
    ],[
        AC_MSG_RESULT(unknown)
        echo "WARNING: You appear to be cross compiling, so there is no way to determine"
        echo "whether your MPFR version is new enough. I am assuming it is."
    ])
    ifelse([$2], , :, [$2])
elif test -n "$mpfr_problem"; then
    AC_MSG_RESULT(problem)
    echo "Sorry, your MPFR version is too old. Disabling."
    ifelse([$3], , :, [$3])
elif test "x$mpfr_found" = "xno" ; then
    AC_MSG_RESULT(not found)
    ifelse([$3], , :, [$3])
fi

AM_CONDITIONAL(LINBOX_HAVE_MPFR, test "x$HAVE_MPFR" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
