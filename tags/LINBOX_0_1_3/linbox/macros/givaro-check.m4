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
dnl Test for Givaro and define GIVARO_CFLAGS and GIVARO_LIBS

AC_DEFUN([LB_CHECK_GIVARO],
[

AC_ARG_WITH(givaro-prefix,[  --with-givaro-prefix=PFX Prefix where GIVARO is installed (optional)],
[givaro_prefix="$withval"],[givaro_prefix=""])

min_givaro_version=ifelse([$1], ,3.0,$1)
AC_MSG_CHECKING(for GIVARO >= $min_givaro_version)

if test x$givaro_prefix = x; then
	givaro_prefix=/usr
else 
	LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${givaro_prefix}/lib" 
	export LD_LIBRARY_PATH
fi

dnl Check for existence

if test "x${givaro_prefix}" != "x/usr" -a "x${givaro_prefix}" != "x/usr/local"; then
	GIVARO_CFLAGS="-I${givaro_prefix}/include"
	GIVARO_LIBS="-L${givaro_prefix}/lib -lgivaro"
else
	GIVARO_CFLAGS=
	GIVARO_LIBS=-lgivaro
fi

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${GIVARO_CFLAGS} ${GMP_CFLAGS}" 
LIBS="${GIVARO_LIBS} ${GMP_LIBS}"

AC_TRY_LINK(
[#include <givinteger.h>],
[Integer a;],
[
AC_TRY_RUN(
[#include <givconfig.h>
int main () {  if (GIVARO_VERSION < 3) return -1; else return 0; }
],[
AC_MSG_RESULT(found)
AC_SUBST(GIVARO_CFLAGS)
AC_SUBST(GIVARO_LIBS)
AC_DEFINE(HAVE_GIVARO,1,[Define if GIVARO is installed])

# N.B. Put real definitions here when we add header files

HAVE_GIVARO=yes

ifelse([$2], , :, [$2])
],[
AC_MSG_RESULT(not found)
echo "Sorry, your GIVARO version is too old. Disabling."

unset GIVARO_CFLAGS
unset GIVARO_LIBS

ifelse([$3], , :, [$3])
],[
AC_MSG_RESULT(unknown)
echo "WARNING: You appear to be cross compiling, so there is no way to determine"
echo "whether your GIVARO version is new enough. I am assuming it is."

HAVE_GIVARO=yes

AC_SUBST(GIVARO_CFLAGS)
AC_SUBST(GIVARO_LIBS)
AC_DEFINE(HAVE_GIVARO,1,[Define if GIVARO is installed])

ifelse([$2], , :, [$2])
])
],
[
AC_MSG_RESULT(not found)

unset GIVARO_CFLAGS
unset GIVARO_LIBS

ifelse([$3], , :, [$3])
])

AM_CONDITIONAL(HAVE_GIVARO, test "x$HAVE_GIVARO" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
