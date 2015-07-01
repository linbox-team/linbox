# Check for Givaro
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

AC_ARG_WITH(givaro-prefix,[  --with-givaro-prefix=PFX      Prefix where GIVARO is installed (optional)],
[givaro_prefix="$withval"],[givaro_prefix=""])

min_givaro_version=ifelse([$1], ,4.0,$1)
AC_MSG_CHECKING(for GIVARO >= $min_givaro_version)

if test x$givaro_prefix != x; then
	export LD_LIBRARY_PATH=$givaro_prefix/lib:$LD_LIBRARY_PATH
	export CPLUS_INCLUDE_PATH=$givaro_prefix/include:$CPLUS_INCLUDE_PATH
fi

dnl Check for existence

AC_CHECK_LIB(givaro, GetTime,
[
dnl Check if the version is new enough
dnl FIXME

GIVARO_CFLAGS="-I$(givaro_prefix)/include"
GIVARO_LIBS="-L$(givaro_prefix)/lib -lgivaro"
AC_SUBST(GIVARO_CFLAGS)
AC_SUBST(GIVARO_LIBS)
AC_DEFINE(HAVE_GIVARO)
AC_MSG_RESULT(found)
ifelse([$2], , :, [$2])
],
[
AC_MSG_RESULT(not found)
ifelse([$3], , :, [$3])
])

])
