# Check for GMP


dnl LB_CHECK_GMP ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for the GNU Multiprecision library and define GMP_CFLAGS and GMP_LIBS

AC_DEFUN([LB_CHECK_GMP],
[

AC_ARG_WITH(gmp-prefix,[  --with-gmp-prefix=PFX   Prefix where GMP is installed],
[gmp_prefix="$withval"],[gmp_prefix=""])

min_gmp_version=ifelse([$1], ,3.1.1,$1)
AC_MSG_CHECKING(for GMP >= $min_gmp_version)

if test x$gmp_prefix = x; then
	gmp_prefix=/usr
else 
	LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${gmp_prefix}/lib"
	export LD_LIBRARY_PATH
fi

dnl Check for existence

if test "x${gmp_prefix}" != "x/usr" -a "x${gmp_prefix}" != "x/usr/local"; then
	if test -d ${gmp_prefix}/include; then
		GMP_CFLAGS="-I${gmp_prefix}/include"
	else
		GMP_CFLAGS="-I${gmp_prefix}"
	fi

	if test -d ${gmp_prefix}/lib; then
		GMP_LIBS="-L${gmp_prefix}/lib -lgmp"
	else
		GMP_LIBS="-L${gmp_prefix} -lgmp"
	fi
else
	GMP_CFLAGS=
	GMP_LIBS=-lgmp
fi

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${CXXFLAGS} ${GMP_CFLAGS}"
LIBS="${LIBS} ${GMP_LIBS}"

AC_TRY_LINK(
[#include <gmp.h>],
[mpz_t a; mpz_init (a);],
[
AC_TRY_RUN(
[#include <gmp.h>
int main () {  if (__GNU_MP_VERSION < 3) return -1; else return 0; }
],[
AC_MSG_RESULT(found)
AC_SUBST(GMP_CFLAGS)
AC_SUBST(GMP_LIBS)
AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])

# See if we are running GMP 4.0
AC_MSG_CHECKING(whether GMP is 4.0 or greater)
AC_TRY_RUN(
[#include <gmp.h>
int main () { if (__GNU_MP_VERSION < 4) return -1; else return 0; }
],[
AC_MSG_RESULT(yes)


# See if GMP was compiled with --enable-cxx
AC_MSG_CHECKING(whether GMP was compiled with --enable-cxx)
AC_TRY_RUN(
[#include <gmpxx.h>
int main () { mpz_class a(2),b(3),c(5); if ( a+b == c ) return 0; else return -1; }
],[
AC_MSG_RESULT(yes)
GMP_VERSION=""
AC_SUBST(GMP_VERSION)
],[
AC_MSG_RESULT(no)
AC_DEFINE(GMP_NO_CXX,1,[Define if GMP has no <gmpxx.h>])
GMP_VERSION="-DGMP_NO_CXX"
AC_SUBST(GMP_VERSION)
],[
dnl This should never happen
AC_MSG_RESULT(no)
])

],[
AC_MSG_RESULT(no)
AC_DEFINE(GMP_VERSION_3,1,[Define if GMP is version 3.xxx])
GMP_VERSION="-DGMP_VERSION_3"
AC_SUBST(GMP_VERSION)
],[
dnl This should never happen
AC_MSG_RESULT(no)
])

ifelse([$2], , :, [$2])
],[
AC_MSG_RESULT(not found)
echo "Sorry, your GMP version is too old. Disabling."

unset GMP_CFLAGS
unset GMP_LIBS

ifelse([$3], , :, [$3])
],[
AC_MSG_RESULT(unknown)
echo "WARNING: You appear to be cross compiling, so there is no way to determine"
echo "whether your GMP version is new enough. I am assuming it is."
AC_SUBST(GMP_CFLAGS)
AC_SUBST(GMP_LIBS)
AC_DEFINE(HAVE_GMP,1,[Define if GMP is installed])

ifelse([$2], , :, [$2])
])
],
[
AC_MSG_RESULT(not found)

unset GMP_CFLAGS
unset GMP_LIBS

ifelse([$3], , :, [$3])
])

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
