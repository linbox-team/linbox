# Check for ATLAS
# Pascal Giorgi , 2003-03-04
# Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey       99-2-5
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor

dnl LB_CHECK_ATLAS ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for ATLAS and define ATLAS_CFLAGS and ATLAS_LIBS

AC_DEFUN([LB_CHECK_ATLAS],
[

AC_ARG_WITH(atlas-prefix,[  --with-atlas-prefix=PFX   Prefix where ATLAS is installed (optional)],
[atlas_prefix="$withval"],[atlas_prefix=""])
AC_ARG_WITH(atlas-arch,[  --with-atlas-arch=ARCH   Architecture where ATLAS has been compiled (optional)],
[atlas_arch="$withval"],[atlas_arch=""])

min_atlas_version=ifelse([$1], ,3.0,$1)
AC_MSG_CHECKING(for ATLAS >= $min_atlas_version)

if test x$atlas_prefix = x; then
	atlas_prefix=/usr
else 
	LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${atlas_prefix}/lib/${atlas_arch}" 
	export LD_LIBRARY_PATH
fi

dnl Check for existence

if test "x${atlas_prefix}" != "x/usr" -a "x${atlas_prefix}" != "x/usr/local"; then
	ATLAS_CFLAGS="-I${atlas_prefix}/include -I${atlas_prefix}/include/${atlas_arch}"
	ATLAS_LIBS="-L${atlas_prefix}/lib/${atlas_arch} -lcblas -latlas"
else
	ATLAS_CFLAGS=
	ATLAS_LIBS="-lcblas -latlas"
fi

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${CXXFLAGS} ${ATLAS_CFLAGS}" 
LIBS="${LIBS} ${ATLAS_LIBS} "

AC_TRY_LINK(
[#include <cblas.h>],
[double a;],
[
AC_TRY_RUN(
[#include <atlas_buildinfo.h>
int main () {  if  (ATL_VERS == "") return -1; else return 0; }
],[
AC_MSG_RESULT(found)
AC_SUBST(ATLAS_CFLAGS)
AC_SUBST(ATLAS_LIBS)
AC_DEFINE(HAVE_ATLAS,1,[Define if ATLAS is installed])
AC_DEFINE(BLAS_AVAILABLE,,[Define if BLAS routines are available])
# N.B. Put real definitions here when we add header files

HAVE_ATLAS=yes


ifelse([$2], , :, [$2])
],[
AC_MSG_RESULT(not found)
echo "Sorry, ATLAS seems not working."

unset ATLAS_CFLAGS
unset ATLAS_LIBS

ifelse([$3], , :, [$3])
],[
AC_MSG_RESULT(unknown)
echo "WARNING: You appear to be cross compiling, so there is no way to determine"
echo "whether your ATLAS version is new enough. I am assuming it is."

HAVE_ATLAS=yes

AC_SUBST(ATLAS_CFLAGS)
AC_SUBST(ATLAS_LIBS)
AC_DEFINE(HAVE_ATLAS,1,[Define if ATLAS is installed])

ifelse([$2], , :, [$2])
])
],
[
AC_MSG_RESULT(not found)

unset ATLAS_CFLAGS
unset ATLAS_LIBS

ifelse([$3], , :, [$3])
])

AM_CONDITIONAL(HAVE_ATLAS, test "x$HAVE_ATLAS" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
