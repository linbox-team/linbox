# Check for SACLIB
# Bradford Hovinen, 2001-06-13
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_SACLIB ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define SACLIB_CFLAGS and SACLIB_LIBS

AC_DEFUN([LB_CHECK_SACLIB],
[

AC_ARG_WITH(saclib-prefix,[  --with-saclib-prefix=PFX      Prefix where SACLIB is installed (optional)],
[saclib_prefix="$withval"],[saclib_prefix=""])

min_saclib_version=ifelse([$1], ,3.1.1,$1)
AC_MSG_CHECKING(for SACLIB >= $min_saclib_version)

if test x$saclib_prefix != x; then
	export LD_LIBRARY_PATH=$saclib_prefix/lib:$LD_LIBRARY_PATH
	export CPLUS_INCLUDE_PATH=$saclib_prefix/include:$CPLUS_INCLUDE_PATH
else
	saclib_prefix=/usr
fi

dnl Check for existence

SACLIB_CFLAGS="-I${saclib_prefix}/include"
SACLIB_LIBS="-L${saclib_prefix}/lib -lsaclib"

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS=${SACLIB_CFLAGS}
LIBS=${SACLIB_LIBS}

AC_TRY_LINK(
[#include <saclib.h>],
[mpz_t a; mpz_init (a);],
[
AC_MSG_RESULT(found)
AC_MSG_CHECKING(SACLIB major version number)
AC_TRY_RUN(
[#include <saclib.h>
int main () {  if (__GNU_MP_VERSION < 3) return -1; else return 0; }
],[
AC_MSG_RESULT(>= 3.0)
AC_SUBST(SACLIB_CFLAGS)
AC_SUBST(SACLIB_LIBS)
AC_DEFINE(HAVE_SACLIB)
],[
AC_MSG_RESULT(< 3.0)
echo "Sorry, your NTL version is too old. Disabling."
])

ifelse([$2], , :, [$2])
],
[
AC_MSG_RESULT(not found)
ifelse([$3], , :, [$3])
])

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}

])
