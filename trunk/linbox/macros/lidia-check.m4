# Check for LIDIA
# Pascal Giorgi, 2001-12-10
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_LIDIA ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for LIDIA Library and define LIDIA_CFLAGS and LIDIA_LIBS

AC_DEFUN([LB_CHECK_LIDIA],
[

AC_ARG_WITH(lidia-prefix,[  --with-lidia-prefix=PFX Prefix where LIDIA is installed (optional)],
[lidia_prefix="$withval"],[lidia_prefix=""])

min_lidia_version=ifelse([$1], ,2.1,$1)
AC_MSG_CHECKING(for LIDIA >= $min_lidia_version)

if test x$lidia_prefix = x; then
	lidia_prefix=/usr/local
else 
	LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${lidia_prefix}/lib"
	export LD_LIBRARY_PATH
fi

dnl Check for existence

if test "x${lidia_prefix}" != "x/usr" -a "x${lidia_prefix}" != "x/usr/local"; then
	LIDIA_CFLAGS="-I${lidia_prefix}/include"
	LIDIA_LIBS="-L${lidia_prefix}/lib -lLiDIA"
else
	LIDIA_CFLAGS=
	LIDIA_LIBS=-lLiDIA
fi

# By default, these should be empty. We set them to include real data
# only if LIDIA is actually found.

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${CXXFLAGS} ${LIDIA_CFLAGS} ${GMP_CFLAGS}"
LIBS="${CXXFLAGS} ${LIDIA_LIBS} ${GMP_LIBS}"

AC_TRY_LINK(
[#include <LiDIA/bigint.h>],
[LiDIA::bigint a;],
[
AC_TRY_RUN(
[#include <LiDIA/LiDIA.h>
#include <iostream>
int main () { if (LIDIA_MAJOR_VERSION < 2) return -1; else return 0; }
],[
AC_MSG_RESULT(found)
AC_SUBST(LIDIA_CFLAGS)
AC_SUBST(LIDIA_LIBS)

AC_DEFINE(HAVE_LIDIA,1,[Define if LiDIA is installed])

# LIDIA was found, so make sure tests and headers get included.

HAVE_LIDIA=yes

ifelse([$2], , :, [$2])
],[
AC_MSG_RESULT(not found)
echo "Sorry, your LIDIA version is too old. Disabling."

unset LIDIA_CFLAGS
unset LIDIA_LIBS

ifelse([$3], , :, [$3])
],[
AC_MSG_RESULT(unknown)
echo "WARNING: You appear to be cross compiling, so there is no way to determine"
echo "whether your LIDIA version is new enough. I am assuming it is."

HAVE_LIDIA=yes

AC_SUBST(LIDIA_CFLAGS)
AC_SUBST(LIDIA_LIBS)
AC_DEFINE(HAVE_LIDIA,1,[Define if LiDIA is installed])

ifelse([$2], , :, [$2])
])
],
[
AC_MSG_RESULT(not found)

unset LIDIA_CFLAGS
unset LIDIA_LIBS

ifelse([$3], , :, [$3])
])

AM_CONDITIONAL(HAVE_LIDIA, test "x$HAVE_LIDIA" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}

])
