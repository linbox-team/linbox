# Check for NTL
# Bradford Hovinen, 2001-06-13
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_NTL ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for Victor Shoup's NTL (Number Theory Library) and define
dnl NTL_CFLAGS and NTL_LIBS

AC_DEFUN([LB_CHECK_NTL],
[

AC_ARG_WITH(ntl-prefix,[  --with-ntl-prefix=PFX      Prefix where NTL is installed (optional)],
[ntl_prefix="$withval"],[ntl_prefix=""])

min_ntl_version=ifelse([$1], ,4.0,$1)
AC_MSG_CHECKING(for NTL >= $min_ntl_version)

if test x$ntl_prefix = x; then
	ntl_prefix=/usr
fi

dnl Check for existence

NTL_CFLAGS="-I${ntl_prefix}/include"
NTL_LIBS="${ntl_prefix}/src/ntl.a"

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS=${NTL_CFLAGS}
LIBS=${NTL_LIBS}

AC_TRY_LINK(
[#include <NTL/ZZ.h>],
[ZZ a;],
[
AC_TRY_RUN(
[#include <NTL/version.h>
#include <iostream>
int main () { if (NTL_MAJOR_VERSION < 4) return -1; else return 0; }
],[
AC_MSG_RESULT(found)
AC_SUBST(NTL_CFLAGS)
AC_SUBST(NTL_LIBS)
AC_DEFINE(HAVE_NTL)

ifelse([$2], , :, [$2])
],[
AC_MSG_RESULT(not found)
echo "Sorry, your NTL version is too old. Disabling."

unset NTL_CFLAGS
unset NTL_LIBS

ifelse([$3], , :, [$3])
])
],
[
AC_MSG_RESULT(not found)
if test x$ntl_prefix != x/usr; then
	AC_MSG_WARN(NTL >= 4.0 was not found. Please double-check the directory you gave.)
fi

unset NTL_CFLAGS
unset NTL_LIBS

ifelse([$3], , :, [$3])
])

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}

])
