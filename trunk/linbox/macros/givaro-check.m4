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
dnl and also for GIVARO_TESTS and GIVARO_HEADERS

AC_DEFUN([LB_CHECK_GIVARO],
[

AC_ARG_WITH(givaro-prefix,[  --with-givaro-prefix=PFX      Prefix where GIVARO is installed (optional)],
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

GIVARO_CFLAGS="-I${givaro_prefix}/include "
GIVARO_LIBS="-L${givaro_prefix}/lib -lgivaro "

# N.B. These should always be empty
GIVARO_TESTS=
GIVARO_HEADERS_BASE=
GIVARO_HEADERS_FIELD= 
GIVARO_HEADERS_BLACKBOX=

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
AC_DEFINE(HAVE_GIVARO)

# N.B. Put real definitions here when we add header files

GIVARO_TESTS="test-givaro-zpz"
GIVARO_HEADERS_BASE="givaro.h"
GIVARO_HEADERS_FIELD="givaro.h givaro-zpz.h givaro-gfq.h"
GIVARO_HEADERS_BLACKBOX=""

ifelse([$2], , :, [$2])
],[
AC_MSG_RESULT(not found)
echo "Sorry, your GIVARO version is too old. Disabling."

unset GIVARO_CFLAGS
unset GIVARO_LIBS


ifelse([$3], , :, [$3])
])
],
[
AC_MSG_RESULT(not found)

unset GIVARO_CFLAGS
unset GIVARO_LIBS

ifelse([$3], , :, [$3])
])

AC_SUBST(GIVARO_TESTS)
AC_SUBST(GIVARO_HEADERS_BASE)
AC_SUBST(GIVARO_HEADERS_FIELD)
AC_SUBST(GIVARO_HEADERS_BLACKBOX)

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
