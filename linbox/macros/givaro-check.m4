# Check for GIVARO
# Bradford Hovinen, 2001-06-13
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_GIVARO ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define GIVARO_CFLAGS and GIVARO_LIBS

AC_DEFUN([LB_CHECK_GIVARO],
[

AC_ARG_WITH(givaro-prefix,[  --with-givaro-prefix=PFX      Prefix where GIVARO is installed (optional)],
[givaro_prefix="$withval"],[givaro_prefix=""])

min_givaro_version=ifelse([$1], ,3.1.1,$1)
AC_MSG_CHECKING(for Givaro >= $min_givaro_version)

if test x$givaro_prefix != x; then
	export LD_LIBRARY_PATH=$givaro_prefix/lib:$LD_LIBRARY_PATH
	export CPLUS_INCLUDE_PATH=$givaro_prefix/include:$CPLUS_INCLUDE_PATH
else
	givaro_prefix=/usr
fi

dnl Check for existence

GIVARO_CFLAGS="-I${givaro_prefix}/include"
GIVARO_LIBS="-L${givaro_prefix}/lib -lgivaro"

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS=${GIVARO_CFLAGS}
LIBS=${GIVARO_LIBS}

AC_TRY_LINK(
[#include <givaro.h>],
[Integer a;],
[
AC_MSG_RESULT(found)
AC_MSG_CHECKING(GIVARO major version number)
AC_TRY_RUN(
[#include <givaro.h>
int main () {  if (GIVARO_VERSION < 3) return -1; else return 0; }
],[
AC_MSG_RESULT(>= 3.0)
AC_SUBST(GIVARO_CFLAGS)
AC_SUBST(GIVARO_LIBS)
AC_DEFINE(HAVE_GIVARO)
],[
AC_MSG_RESULT(< 3.0)
echo "Sorry, your Givaro version is too old. Disabling."
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
